"""Read alignment implementations and alignment parsing."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

from .models import Alignment
from .sequence_io import read_fastq, reverse_complement, write_reference_fasta


def parse_cigar(cigar: str) -> List[Tuple[str, int]]:
    """Parse a SAM CIGAR string into (operation, length) tuples."""
    parts: List[Tuple[str, int]] = []
    number = ""
    for char in cigar:
        if char.isdigit():
            number += char
        else:
            if not number:
                raise ValueError(f"Invalid CIGAR string: {cigar!r}")
            parts.append((char, int(number)))
            number = ""
    if number:
        raise ValueError(f"Invalid CIGAR string: {cigar!r}")
    return parts


def align_reads_with_minimap2(
    reads_path: str,
    reference_name: str,
    reference_seq: str,
    threads: int,
) -> Dict[str, Alignment]:
    """Run minimap2 and return the best primary alignment per read."""
    if shutil.which("minimap2") is None:
        raise SystemExit("minimap2 was not found on PATH. Install it or use --aligner smithwaterman.")

    best_by_read: Dict[str, Alignment] = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_fasta = f"{tmpdir}/reference.fasta"
        sam_path = f"{tmpdir}/alignments.sam"
        write_reference_fasta(ref_fasta, reference_name, reference_seq)

        cmd = [
            "minimap2",
            "-ax",
            "map-ont",
            "--eqx",
            "--secondary=no",
            "-t",
            str(max(1, threads)),
            ref_fasta,
            reads_path,
        ]
        with open(sam_path, "w") as sam_out:
            subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, text=True, check=True)

        with open(sam_path) as handle:
            for raw_line in handle:
                if not raw_line or raw_line.startswith("@"):
                    continue
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    continue

                read_id = fields[0]
                flag = int(fields[1])
                ref_name_field = fields[2]
                if ref_name_field == "*":
                    continue
                ref_pos_1based = int(fields[3])
                cigar = fields[5]
                seq = fields[9]
                if cigar == "*" or seq == "*":
                    continue
                if flag & 0x900:
                    continue

                orientation = "-" if (flag & 16) else "+"
                ref_pos = ref_pos_1based - 1
                query_pos = 0
                matches = 0
                mismatches = 0
                insertions = 0
                deletions = 0
                aligned_pairs: List[Tuple[Optional[int], Optional[int], str]] = []
                ref_start = ref_pos

                for op, length in parse_cigar(cigar):
                    if op in {"=", "X", "M"}:
                        for _ in range(length):
                            if op == "=":
                                matches += 1
                                pair_op = "M"
                            elif op == "X":
                                mismatches += 1
                                pair_op = "X"
                            else:
                                matches += 1
                                pair_op = "M"
                            aligned_pairs.append((ref_pos, query_pos, pair_op))
                            ref_pos += 1
                            query_pos += 1
                    elif op == "I":
                        for _ in range(length):
                            insertions += 1
                            aligned_pairs.append((None, query_pos, "I"))
                            query_pos += 1
                    elif op in {"D", "N"}:
                        for _ in range(length):
                            deletions += 1
                            aligned_pairs.append((ref_pos, None, "D"))
                            ref_pos += 1
                    elif op == "S":
                        query_pos += length
                    elif op in {"H", "P"}:
                        continue
                    else:
                        raise ValueError(f"Unsupported CIGAR op {op!r} in {cigar!r}")

                ref_end = ref_pos
                score = 0
                for tag in fields[11:]:
                    if tag.startswith("AS:i:"):
                        score = int(tag.split(":", 2)[2])
                        break
                if score == 0:
                    score = matches

                aln = Alignment(
                    read_id=read_id,
                    orientation=orientation,
                    score=score,
                    ref_start=ref_start,
                    ref_end=ref_end,
                    matches=matches,
                    mismatches=mismatches,
                    insertions=insertions,
                    deletions=deletions,
                    aligned_pairs=aligned_pairs,
                )
                existing = best_by_read.get(read_id)
                if existing is None or aln.score > existing.score:
                    best_by_read[read_id] = aln
    return best_by_read


def smith_waterman(
    read_id: str,
    read_seq: str,
    ref_seq: str,
    orientation: str,
    match_score: int,
    mismatch_penalty: int,
    gap_penalty: int,
) -> Alignment:
    """Local alignment of read_seq to ref_seq using Smith-Waterman."""
    n = len(read_seq)
    m = len(ref_seq)
    prev = [0] * (m + 1)
    curr = [0] * (m + 1)
    traceback: List[bytearray] = [bytearray(m + 1) for _ in range(n + 1)]

    best_score = 0
    best_i = 0
    best_j = 0

    for i in range(1, n + 1):
        curr[0] = 0
        read_base = read_seq[i - 1]
        for j in range(1, m + 1):
            ref_base = ref_seq[j - 1]
            diag_score = prev[j - 1] + (match_score if read_base == ref_base else -mismatch_penalty)
            up_score = prev[j] - gap_penalty
            left_score = curr[j - 1] - gap_penalty
            score = max(0, diag_score, up_score, left_score)
            curr[j] = score
            if score == 0:
                traceback[i][j] = 0
            elif score == diag_score:
                traceback[i][j] = 1  # diagonal
            elif score == up_score:
                traceback[i][j] = 2  # insertion relative to reference
            else:
                traceback[i][j] = 3  # deletion relative to reference
            if score > best_score:
                best_score = score
                best_i = i
                best_j = j
        prev, curr = curr, prev

    i = best_i
    j = best_j
    ref_end = best_j
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0
    pairs_reversed: List[Tuple[Optional[int], Optional[int], str]] = []

    while i > 0 and j > 0 and traceback[i][j] != 0:
        direction = traceback[i][j]
        if direction == 1:
            ref_pos = j - 1
            read_pos = i - 1
            if read_seq[read_pos] == ref_seq[ref_pos]:
                matches += 1
                op = "M"
            else:
                mismatches += 1
                op = "X"
            pairs_reversed.append((ref_pos, read_pos, op))
            i -= 1
            j -= 1
        elif direction == 2:
            insertions += 1
            pairs_reversed.append((None, i - 1, "I"))
            i -= 1
        else:
            deletions += 1
            pairs_reversed.append((j - 1, None, "D"))
            j -= 1

    pairs_reversed.reverse()
    return Alignment(
        read_id=read_id,
        orientation=orientation,
        score=best_score,
        ref_start=j,
        ref_end=ref_end,
        matches=matches,
        mismatches=mismatches,
        insertions=insertions,
        deletions=deletions,
        aligned_pairs=pairs_reversed,
    )


def best_alignment_for_read(
    read_id: str,
    read_seq: str,
    ref_seq: str,
    match_score: int,
    mismatch_penalty: int,
    gap_penalty: int,
) -> Alignment:
    forward = smith_waterman(
        read_id, read_seq, ref_seq, "+", match_score, mismatch_penalty, gap_penalty
    )
    reverse = smith_waterman(
        read_id, reverse_complement(read_seq), ref_seq, "-", match_score, mismatch_penalty, gap_penalty
    )
    return forward if forward.score >= reverse.score else reverse


def collect_alignments(
    reads_path: str,
    ref_name: str,
    ref_seq: str,
    aligner: str,
    threads: int,
    match_score: int,
    mismatch_penalty: int,
    gap_penalty: int,
) -> List[Alignment]:
    if aligner == "minimap2":
        best_by_read = align_reads_with_minimap2(
            reads_path=reads_path,
            reference_name=ref_name,
            reference_seq=ref_seq,
            threads=threads,
        )
        return list(best_by_read.values())

    alignments: List[Alignment] = []
    for read_id, read_seq in read_fastq(reads_path):
        if not read_seq:
            continue
        aln = best_alignment_for_read(
            read_id=read_id,
            read_seq=read_seq,
            ref_seq=ref_seq,
            match_score=match_score,
            mismatch_penalty=mismatch_penalty,
            gap_penalty=gap_penalty,
        )
        alignments.append(aln)
    return alignments
