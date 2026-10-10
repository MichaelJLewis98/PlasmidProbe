"""Calculate feature coverage and evidence from read alignments."""

from __future__ import annotations

from typing import Dict, List

from .alignments import collect_alignments
from .models import Feature, FeatureResult


def analyse_features(
    reads_path: str,
    ref_name: str,
    ref_seq: str,
    features: List[Feature],
    aligner: str,
    threads: int,
    match_score: int,
    mismatch_penalty: int,
    gap_penalty: int,
    min_alignment_score: int,
    min_alignment_score_density: float,
    min_read_identity: float,
    min_feature_breadth: float,
    min_feature_depth: float,
    min_evidence_reads: int,
) -> List[FeatureResult]:
    ref_len = len(ref_seq)
    depth = [0] * ref_len
    identity_sum = [0.0] * ref_len
    identity_count = [0] * ref_len
    feature_evidence_reads: Dict[str, set[str]] = {feature.name: set() for feature in features}

    alignments = collect_alignments(
        reads_path=reads_path,
        ref_name=ref_name,
        ref_seq=ref_seq,
        aligner=aligner,
        threads=threads,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_penalty=gap_penalty,
    )

    for aln in alignments:
        if (
            aln.score < min_alignment_score
            or aln.score_density < min_alignment_score_density
            or aln.identity < min_read_identity
        ):
            continue

        read_covered_positions = set()
        read_matching_positions = set()
        for ref_pos, _read_pos, op in aln.aligned_pairs:
            if ref_pos is None:
                continue
            if 0 <= ref_pos < ref_len:
                depth[ref_pos] += 1
                identity_count[ref_pos] += 1
                if op == "M":
                    identity_sum[ref_pos] += 1.0
                    read_matching_positions.add(ref_pos)
                read_covered_positions.add(ref_pos)

        for feature in features:
            covered_in_feature = [pos for pos in read_covered_positions if feature.start <= pos < feature.end]
            if not covered_in_feature:
                continue
            matching_in_feature = [pos for pos in covered_in_feature if pos in read_matching_positions]
            read_feature_breadth = len(set(covered_in_feature)) / (feature.end - feature.start)
            read_feature_identity = len(matching_in_feature) / len(set(covered_in_feature))
            if read_feature_breadth >= min_feature_breadth and read_feature_identity >= min_read_identity:
                feature_evidence_reads[feature.name].add(aln.read_id)

    results: List[FeatureResult] = []
    for feature in features:
        feature_depths = depth[feature.start : feature.end]
        covered_bases = sum(1 for d in feature_depths if d > 0)
        length = feature.end - feature.start
        breadth = covered_bases / length if length else 0.0
        mean_depth = sum(feature_depths) / length if length else 0.0

        per_base_identities = []
        for pos in range(feature.start, feature.end):
            if identity_count[pos] > 0:
                per_base_identities.append(identity_sum[pos] / identity_count[pos])
        mean_identity = sum(per_base_identities) / len(per_base_identities) if per_base_identities else 0.0

        evidence_reads = len(feature_evidence_reads[feature.name])
        present = (
            breadth >= min_feature_breadth
            and mean_depth >= min_feature_depth
            and evidence_reads >= min_evidence_reads
        )
        results.append(
            FeatureResult(
                feature=feature.name,
                source_name=feature.source_name or feature.name,
                start=feature.start + 1,
                end=feature.end,
                length=length,
                status="EVIDENCE_DETECTED" if present else "NOT_DETECTED",
                covered_bases=covered_bases,
                breadth=breadth,
                mean_depth=mean_depth,
                mean_identity=mean_identity,
                evidence_reads=evidence_reads,
            )
        )
    return results
