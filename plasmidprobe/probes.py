"""Generate overlapping probes and check their k-mers in reads."""

from __future__ import annotations

from typing import Dict, List

from .models import Feature, Probe
from .sequence_io import read_fastq, reverse_complement


def generate_probes(
    reference_seq: str,
    features: List[Feature],
    probe_length: int = 50,
    probe_step: int = 25,
    ) -> List[Probe]:
    probes: List[Probe] = []
    for feature in features:
        feature_seq = reference_seq[feature.start:feature.end]

        if len(feature_seq) < probe_length:
            probes.append(
                Probe(
                    feature_name=feature.name,
                    source_name=feature.source_name or feature.name,
                    start=feature.start,
                    end=feature.end,
                    sequence=feature_seq,
                )
            )
            continue

        for i in range(
            0,
            len(feature_seq) - probe_length + 1,
            probe_step,
        ):
            probe_seq = feature_seq[i : i + probe_length]

            probes.append(
                Probe(
                    feature_name=feature.name,
                    source_name=feature.source_name or feature.name,
                    start=feature.start + i,
                    end=feature.start + i + probe_length,
                    sequence=probe_seq,
                )
            )

    return probes


def probe_hits_from_reads(
    reads_path: str,
    probes: List[Probe],
    k: int = 15,
    ) -> Dict[str, int]:
    def kmers(seq: str, k: int):
        return {
            seq[i:i+k]
            for i in range(len(seq) - k + 1)
        }

    read_kmers = set()

    for _, read_seq in read_fastq(reads_path):
        if len(read_seq) >= k:
            read_kmers.update(kmers(read_seq, k))
            read_kmers.update(
                kmers(reverse_complement(read_seq), k)
            )

    hits: Dict[str, int] = {}

    for probe in probes:

        if len(probe.sequence) < k:
            continue

        probe_kmers = kmers(probe.sequence, k)

        overlap = len(
            probe_kmers.intersection(read_kmers)
        )

        if overlap >= max(1, int(0.3 * len(probe_kmers))):
            hits[probe.feature_name] = (
                hits.get(probe.feature_name, 0) + 1
            )

    return hits
