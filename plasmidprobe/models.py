"""Data structures shared by the analysis modules."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass
class Feature:
    name: str
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive
    strand: str = "+"
    source_name: str = ""


@dataclass
class Alignment:
    read_id: str
    orientation: str
    score: int
    ref_start: int
    ref_end: int
    matches: int
    mismatches: int
    insertions: int
    deletions: int
    aligned_pairs: List[Tuple[Optional[int], Optional[int], str]]

    @property
    def aligned_length(self) -> int:
        return self.matches + self.mismatches + self.insertions + self.deletions

    @property
    def identity(self) -> float:
        denom = self.matches + self.mismatches + self.insertions + self.deletions
        return self.matches / denom if denom else 0.0

    @property
    def score_density(self) -> float:
        return self.score / self.aligned_length if self.aligned_length else 0.0


@dataclass
class FeatureResult:
    feature: str
    source_name: str
    start: int
    end: int
    length: int
    status: str
    covered_bases: int
    breadth: float
    mean_depth: float
    mean_identity: float
    evidence_reads: int


@dataclass
class Probe:
    feature_name: str
    source_name: str
    start: int
    end: int
    sequence: str
