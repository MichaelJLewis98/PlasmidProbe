"""Read reference, read, and feature annotation files."""

from __future__ import annotations

import csv
import gzip
import re
from typing import Dict, Iterator, List, Optional, Tuple

from Bio import SeqIO

from .models import Feature

DNA_RE = re.compile(r"[^ACGTNacgtn]")
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def clean_seq(seq: str) -> str:
    return DNA_RE.sub("", seq).upper()


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()


def read_fasta(path: str) -> Dict[str, str]:
    records: Dict[str, List[str]] = {}
    name = None
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                if not name:
                    raise ValueError("Found a FASTA header with no name.")
                if name in records:
                    raise ValueError(f"Duplicate FASTA record name: {name}")
                records[name] = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence encountered before first header.")
                records[name].append(clean_seq(line))
    return {name: "".join(parts) for name, parts in records.items()}


def read_single_fasta(path: str) -> Tuple[str, str]:
    records = read_fasta(path)
    if not records:
        raise ValueError("No reference sequence found in FASTA file.")
    if len(records) > 1:
        raise ValueError("Reference FASTA should contain exactly one plasmid sequence.")
    return next(iter(records.items()))


def read_fastq(path: str) -> Iterator[Tuple[str, str]]:
    with open_text(path) as handle:
        while True:
            header = handle.readline().rstrip()
            if not header:
                break
            seq = handle.readline().rstrip()
            plus = handle.readline().rstrip()
            qual = handle.readline().rstrip()
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"FASTQ parsing failed near header: {header}")
            read_id = header[1:].split()[0]
            yield read_id, clean_seq(seq)


def parse_coordinate_values(value_text: str) -> List[int]:
    """Parse either a single coordinate or a brace-delimited list of coordinates."""
    cleaned = value_text.strip().replace("{", "").replace("}", "")
    values = []
    for part in cleaned.split(","):
        part = part.strip()
        if part:
            values.append(int(part))
    if not values:
        raise ValueError(f"Could not parse coordinate value from: {value_text!r}")
    return values


def read_features_csv(path: str, reference_length: int) -> List[Feature]:
    features: List[Feature] = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"name", "start", "end"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Feature CSV is missing required column(s): {', '.join(sorted(missing))}")
        for row in reader:
            name = row["name"].strip()
            start_1based = int(row["start"])
            end_1based = int(row["end"])
            strand = row.get("strand", "+").strip() or "+"
            if not name:
                raise ValueError("Feature with blank name found.")
            if start_1based < 1 or end_1based < start_1based or end_1based > reference_length:
                raise ValueError(
                    f"Invalid coordinates for {name}: start={start_1based}, end={end_1based}, "
                    f"reference length={reference_length}"
                )
            features.append(
                Feature(
                    name=name,
                    start=start_1based - 1,
                    end=end_1based,
                    strand=strand,
                    source_name=name,
                )
            )
    return features


def read_features_index_txt(path: str, reference_length: int) -> List[Feature]:
    """Read the text index format with >feature headers and start/end lines."""
    features: List[Feature] = []
    current_name: Optional[str] = None
    coord_re = re.compile(r"start\s+(.+?),\s*end\s+(.+)$", re.IGNORECASE)

    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].strip()
                if not current_name:
                    raise ValueError("Found a blank feature name in index text file.")
                continue

            if current_name is None:
                raise ValueError(f"Coordinate line appears before a feature name: {line}")

            match = coord_re.match(line)
            if not match:
                raise ValueError(
                    f"Could not parse coordinate line for {current_name!r}: {line!r}. "
                    "Expected format like: start 8, end 883"
                )

            starts = parse_coordinate_values(match.group(1))
            ends = parse_coordinate_values(match.group(2))
            if len(starts) != len(ends):
                raise ValueError(
                    f"Feature {current_name!r} has {len(starts)} starts but {len(ends)} ends."
                )

            for index, (start_1based, end_1based) in enumerate(zip(starts, ends), start=1):
                if start_1based < 1 or end_1based < start_1based or end_1based > reference_length:
                    raise ValueError(
                        f"Invalid coordinates for {current_name}: start={start_1based}, "
                        f"end={end_1based}, reference length={reference_length}"
                    )
                feature_name = current_name if len(starts) == 1 else f"{current_name}_part{index}"
                features.append(
                    Feature(
                        name=feature_name,
                        start=start_1based - 1,
                        end=end_1based,
                        strand="+",
                        source_name=current_name,
                    )
                )
    return features


def _genbank_feature_name(feat) -> str:
    """Pick a human-readable name from common GenBank qualifiers."""
    for key in ("label", "gene", "locus_tag", "product", "note"):
        vals = feat.qualifiers.get(key)
        if vals:
            value = str(vals[0]).strip()
            if value:
                return value
    return feat.type


def read_features_genbank(path: str, reference_length: int) -> List[Feature]:
    """
    Read annotated features from a GenBank file.

    Uses SeqIO.read because your file should contain one plasmid record.
    """
    record = SeqIO.read(path, "genbank")
    if len(record.seq) != reference_length:
        raise ValueError(
            f"GenBank sequence length ({len(record.seq)}) does not match "
            f"reference FASTA length ({reference_length})."
        )

    features: List[Feature] = []
    seen_names: Dict[str, int] = {}

    for feat in record.features:
        if feat.type == "source":
            continue

        name = _genbank_feature_name(feat)
        if not name:
            continue

        seen_names[name] = seen_names.get(name, 0) + 1
        unique_base_name = name if seen_names[name] == 1 else f"{name}_{seen_names[name]}"

        parts = getattr(feat.location, "parts", [feat.location])
        for idx, part in enumerate(parts, start=1):
            start_1based = int(part.start) + 1
            end_1based = int(part.end)
            strand = "-" if part.strand == -1 else "+"

            if start_1based < 1 or end_1based > reference_length:
                raise ValueError(
                    f"Invalid GenBank coordinates for {name}: "
                    f"start={start_1based}, end={end_1based}, "
                    f"reference length={reference_length}"
                )

            feature_name = unique_base_name if len(parts) == 1 else f"{unique_base_name}_part{idx}"
            features.append(
                Feature(
                    name=feature_name,
                    start=start_1based - 1,
                    end=end_1based,
                    strand=strand,
                    source_name=name,
                )
            )

    return features


def read_features(path: str, reference_length: int) -> List[Feature]:
    """Read feature annotations from CSV, index.txt, or GenBank."""
    lower = path.lower()
    if lower.endswith((".gb", ".gbk", ".genbank")):
        return read_features_genbank(path, reference_length)

    with open(path) as handle:
        first_nonempty = ""
        for line in handle:
            if line.strip():
                first_nonempty = line.strip()
                break

    if first_nonempty.startswith(">"):
        return read_features_index_txt(path, reference_length)

    return read_features_csv(path, reference_length)


def write_reference_fasta(path: str, reference_name: str, reference_seq: str) -> None:
    with open(path, "w") as handle:
        handle.write(f">{reference_name}\n")
        for i in range(0, len(reference_seq), 80):
            handle.write(reference_seq[i : i + 80] + "\n")
