"""Command-line interface for PlasmidProbe."""

from __future__ import annotations

import argparse

from .analysis import analyse_features
from .features import extract_target_feature_reference
from .probes import generate_probes, probe_hits_from_reads
from .results import write_results
from .sequence_io import read_features, read_single_fasta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Assess evidence for annotated plasmid features in sequencing reads."
    )
    parser.add_argument("--reads", required=True, help="Input sequencing reads in FASTQ or FASTQ.GZ format.")
    parser.add_argument(
        "--reference",
        required=True,
        help="Single-record reference FASTA; feature mode extracts its target from this sequence.",
    )
    parser.add_argument(
        "--features",
        help="Feature annotations in GenBank, CSV (name,start,end,strand), or text index format.",
    )
    parser.add_argument(
        "--mode",
        choices=["plasmid", "feature", "probes"],
        default="plasmid",
        help="Score all annotated regions (plasmid), one named region (feature), or quick sequence probes (probes).",
    )
    parser.add_argument(
        "--target-name",
        default=None,
        help="Name of the annotated feature to extract in feature mode (required in that mode).",
    )
    parser.add_argument(
        "--aligner",
        choices=["smithwaterman", "minimap2"],
        default="smithwaterman",
        help="Choose smithwaterman for the built-in Python aligner or minimap2 for external long-read mapping.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use with minimap2.",
    )
    parser.add_argument(
        "--probe-length",
        type=int,
        default=50,
        help="Length of each probe in probe mode.",
    )

    parser.add_argument(
        "--probe-step",
        type=int,
        default=25,
        help="Distance between the starts of adjacent probes.",
    )
    parser.add_argument(
        "--out",
        default="results/plasmid_feature_results.csv",
        help="Output CSV file (default: results/plasmid_feature_results.csv).",
    )
    parser.add_argument("--match", type=int, default=2, help="Alignment match score.")
    parser.add_argument("--mismatch", type=int, default=3, help="Alignment mismatch penalty.")
    parser.add_argument("--gap", type=int, default=5, help="Alignment insertion/deletion gap penalty.")
    parser.add_argument("--min-alignment-score", type=int, default=40, help="Discard read alignments below this score.")
    parser.add_argument(
        "--min-alignment-score-density",
        type=float,
        default=0.0,
        help="Discard read alignments whose score per aligned base is below this value.",
    )
    parser.add_argument("--min-read-identity", type=float, default=0.75, help="Discard reads below this alignment identity.")
    parser.add_argument("--min-feature-breadth", type=float, default=0.80, help="Minimum fraction of feature bases covered.")
    parser.add_argument("--min-feature-depth", type=float, default=1.0, help="Minimum mean read depth across a feature.")
    parser.add_argument(
        "--min-evidence-reads",
        "--min-supporting-reads",
        dest="min_evidence_reads",
        type=int,
        default=1,
        help="Minimum reads providing feature-level evidence.",
    )
    args = parser.parse_args()

    ref_name, ref_seq = read_single_fasta(args.reference)
    if args.mode == "plasmid":
        if not args.features:
            raise SystemExit("--features is required in plasmid mode.")
        features = read_features(args.features, len(ref_seq))
        if not features:
            raise SystemExit("No features were found in the feature annotation file.")
    elif args.mode == "probes":
        if not args.features:
            raise SystemExit(
                "--features required in probe mode."
            )

        features = read_features(
            args.features,
            len(ref_seq),
        )

        probes = generate_probes(
            reference_seq=ref_seq,
            features=features,
            probe_length=args.probe_length,
            probe_step=args.probe_step,
        )

        hits = probe_hits_from_reads(
            args.reads,
            probes,
        )

        print()
        print("Probe-based evidence")
        print("-------------------")

        feature_probe_counts = {}

        for probe in probes:
            feature_probe_counts[probe.feature_name] = (
                feature_probe_counts.get(
                    probe.feature_name,
                    0,
                )
                + 1
            )

        for feature_name, total in feature_probe_counts.items():

            detected = hits.get(feature_name, 0)

            fraction = detected / total

            status = (
                "EVIDENCE_DETECTED"
                if fraction >= 0.2
                else "NOT_DETECTED"
            )

            print(
                f"{feature_name}: "
                f"{detected}/{total} probes "
                f"({fraction:.1%}) "
                f"{status}"
            )

        return
    else:
        if not args.features:
            raise SystemExit("--features is required in feature mode so the target sequence can be extracted.")
        all_features = read_features(args.features, len(ref_seq))
        target_feature_name, target_seq, features = extract_target_feature_reference(
            reference_name=ref_name,
            reference_seq=ref_seq,
            features=all_features,
            target_name=args.target_name,
        )
        ref_name = target_feature_name
        ref_seq = target_seq

    results = analyse_features(
        reads_path=args.reads,
        ref_name=ref_name,
        ref_seq=ref_seq,
        features=features,
        aligner=args.aligner,
        threads=args.threads,
        match_score=args.match,
        mismatch_penalty=args.mismatch,
        gap_penalty=args.gap,
        min_alignment_score=args.min_alignment_score,
        min_alignment_score_density=args.min_alignment_score_density,
        min_read_identity=args.min_read_identity,
        min_feature_breadth=args.min_feature_breadth,
        min_feature_depth=args.min_feature_depth,
        min_evidence_reads=args.min_evidence_reads,
    )
    write_results(results, args.out)

    print(f"Reference: {ref_name} ({len(ref_seq)} bp)")
    print(
        f"Mode: {args.mode} | Aligner: {args.aligner} | "
        f"Analysed {len(features)} feature(s). Results written to: {args.out}"
    )
    for result in results:
        print(
            f"{result.feature}: {result.status} | "
            f"coordinates={result.start}-{result.end} | "
            f"breadth={result.breadth:.1%} | "
            f"mean depth={result.mean_depth:.2f} | "
            f"mean identity={result.mean_identity:.1%} | "
            f"evidence reads={result.evidence_reads}"
        )
