"""Operations on annotated plasmid features."""

from __future__ import annotations

from typing import List, Optional, Tuple

from .models import Feature
from .sequence_io import reverse_complement


def extract_target_feature_reference(
    reference_name: str,
    reference_seq: str,
    features: List[Feature],
    target_name: Optional[str],
) -> Tuple[str, str, List[Feature]]:
    """Extract one annotated feature sequence from the reference plasmid."""
    if not target_name:
        raise SystemExit("--target-name is required in feature mode.")

    exact_matches = [feature for feature in features if feature.name == target_name]
    source_matches = [feature for feature in features if feature.source_name == target_name]

    if exact_matches:
        matches = exact_matches
    elif source_matches:
        matches = source_matches
    else:
        raise SystemExit(f"Target feature {target_name!r} was not found in the annotation file.")

    if len(matches) > 1:
        raise SystemExit(
            f"Target feature {target_name!r} resolves to multiple intervals in the annotation file. "
            "Use a single contiguous feature or run in plasmid mode instead."
        )

    target_feature = matches[0]
    target_seq = reference_seq[target_feature.start : target_feature.end]
    if target_feature.strand == "-":
        target_seq = reverse_complement(target_seq)

    target_feature_name = target_feature.name
    target_features = [
        Feature(
            name=target_feature_name,
            start=0,
            end=len(target_seq),
            strand=target_feature.strand,
            source_name=target_feature.source_name or target_feature.name,
        )
    ]
    return target_feature_name, target_seq, target_features
