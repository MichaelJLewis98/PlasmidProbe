"""PlasmidProbe reference-guided plasmid feature analysis."""

from .analysis import analyse_features
from .alignments import best_alignment_for_read, collect_alignments, parse_cigar, smith_waterman
from .features import extract_target_feature_reference
from .models import Alignment, Feature, FeatureResult, Probe
from .probes import generate_probes, probe_hits_from_reads
from .results import write_results
from .sequence_io import (
    clean_seq,
    open_text,
    parse_coordinate_values,
    read_fasta,
    read_fastq,
    read_features,
    read_features_csv,
    read_features_genbank,
    read_features_index_txt,
    read_single_fasta,
    reverse_complement,
)

__all__ = [
    "Alignment", "Feature", "FeatureResult", "Probe", "analyse_features",
    "best_alignment_for_read", "clean_seq", "collect_alignments",
    "extract_target_feature_reference", "generate_probes", "open_text",
    "parse_cigar", "parse_coordinate_values", "probe_hits_from_reads",
    "read_fasta", "read_fastq", "read_features", "read_features_csv",
    "read_features_genbank", "read_features_index_txt", "read_single_fasta",
    "reverse_complement", "smith_waterman", "write_results",
]
