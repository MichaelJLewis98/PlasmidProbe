# PlasmidProbe

**Reference-guided checks for evidence of annotated plasmid features in sequencing reads.**

PlasmidProbe is a Python-based CLI tool for identifying evidence of expected plasmid features in sequencing reads. It compares raw reads with an annotated plasmid reference using alignment according to short, overlapping sequence probes. This is a useful diagnostic to assess evidence for genes, resistance markers, origins of replication, regulatory elements, and any other annotated features you desire.

The project grew out of frustrating DNA Mini-prep samples extracted from Cupriavidus necator cells where the DNA seemed of good quality, but the sequencing provider could not produce a plasmid map because the reads had too little coverage, were too short, or did not meet its quality threshold. It is intended as a secondary check for expected DNA genetic-elements, genes, inserts, and backbone regions. In my case, I was assembling plasmids with many inserts with  modular cloning (MoClo) using GoldenGate reactions. I hope this helps others who receive disapointing sequencing results decide whether to transform new cells or send new samples for sequencing. Missing features after transformation could be due to the host cell's restriction-modification system.

PlasmidProbe reports evidence at the annotated-feature level. A matching fragment does not show that the complete feature is intact or establish the structure of the full plasmid, just that the DNA is there. Its minimap2 option currently uses the `map-ont` preset; the built-in aligner's scoring values can be adjusted. `EVIDENCE_DETECTED` means the configured thresholds were met; `NOT_DETECTED` means they were not according to the threshold, so pick your own carefully.

## What it does

- Reads FASTQ and gzip-compressed FASTQ files.
- Uses a single-record FASTA as the reference sequence.
- Reads feature annotations from GenBank, CSV, or the project's simple text index format.
- Aligns reads in either orientation using a built-in Smith–Waterman implementation or, optionally, minimap2.
- Reports per-feature breadth of coverage, mean depth, mean identity, the number of reads meeting feature-level criteria, and a threshold-based status in CSV format.
- Offers a probe mode that checks for short sequence matches and prints a summary to the terminal.

The built-in aligner is intended for small plasmids and modest read sets. Its runtime can grow substantially with read length, reference length, and read count. Use minimap2 for larger inputs.

## Requirements

- Python 3.9 or newer
- Biopython (installed from `requirements.txt`)
- Optional: [minimap2](https://github.com/lh3/minimap2) installed separately and available on `PATH` when using `--aligner minimap2`

Create an environment and install the Python dependency:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Code layout

The `plasmidprobe/` package separates file handling, alignment, feature scoring, probe searches, results, and the cli. Run it with `python plasmid_probe.py` or `python -m plasmidprobe`.

## Local file layout

Keep local input data and generated outputs in the following folders. They are excluded from Git by `.gitignore`.

```text
data/
├── raw_reads/          # FASTQ/FASTQ.GZ reads
├── references/         # one-record plasmid FASTA and/or GenBank files
└── feature_annotations/ # CSV or text-index annotations
results/                # generated CSV results and intermediate files
```

The existing local files have been grouped by role in these folders. Their original sequence and result filenames have been retained. The feature index files have been given more descriptive names.

## Input formats

The reference FASTA must contain only one sequence. The feature file can be:

- GenBank (`.gb`, `.gbk`, or `.genbank`) with the same sequence length as the reference FASTA.
- CSV with at least `name,start,end` columns and an optional `strand` column. Coordinates are 1-based and inclusive.
- The project's text index format, with a feature name on a `>` line followed by a line such as `start 1050, end 1875`. Coordinates are 1-based and inclusive; brace-delimited coordinate lists are supported for multi-interval features.

## Usage

### Score annotated features across a candidate plasmid

```bash
python plasmid_probe.py \
  --mode plasmid \
  --aligner minimap2 \
  --threads 4 \
  --reads data/raw_reads/sample.fastq.gz \
  --reference data/references/candidate_plasmid.fasta \
  --features data/feature_annotations/candidate_plasmid.gb \
  --out results/sample_feature_evidence.csv
```

For a small input, omit `--aligner minimap2` to use the built-in Smith–Waterman aligner. Defaults include a minimum alignment score of 40, minimum read identity of 0.75, minimum feature breadth of 0.80, minimum mean feature depth of 1, and minimum evidence-read count of 1. Adjust these with `--min-alignment-score`, `--min-read-identity`, `--min-feature-breadth`, `--min-feature-depth`, and `--min-evidence-reads` as appropriate for the experiment. `--min-supporting-reads` remains as an alias for older commands.

### Check one annotated feature

Feature mode extracts one named feature from the annotated reference before aligning reads:

```bash
python plasmid_probe.py \
  --mode feature \
  --aligner minimap2 \
  --threads 4 \
  --reads data/raw_reads/sample.fastq.gz \
  --reference data/references/candidate_plasmid.fasta \
  --features data/feature_annotations/candidate_plasmid.gb \
  --target-name EFE_CDS \
  --out results/sample_efe_evidence.csv
```

### Run a quick probe check

Probe mode creates overlapping probes from the annotated features, checks their 15-mers against the reads in both orientations, and prints the detected fraction for each feature:

```bash
python plasmid_probe.py \
  --mode probes \
  --reads data/raw_reads/sample.fastq.gz \
  --reference data/references/candidate_plasmid.fasta \
  --features data/feature_annotations/candidate_plasmid.gb
```

Probe mode is a lightweight quick sequence-match screen. `EVIDENCE_DETECTED` requires at least 20% of a feature's probes to hit, so treat that result as preliminary evidence.

## Demo and tests

The files in `tests/fixtures/` are short synthetic sequences made for this repository. Give them a try. Three of the four annotated regions have a synthetic read, while one is intentionally missing/uncovered.

Run the demo with:

```bash
python plasmid_probe.py \
  --mode plasmid \
  --reads tests/fixtures/demo_reads.fastq \
  --reference tests/fixtures/demo_reference.fasta \
  --features tests/fixtures/demo_features.csv \
  --out results/demo.csv
```

The output should mark `demo_feature_1`, `demo_feature_2`, and `demo_feature_4` as `EVIDENCE_DETECTED`; `demo_feature_3_unobserved` should be `NOT_DETECTED`.

The tests use Python's built-in `unittest` runner. Install the project requirements first, then run:

```bash
python -m unittest discover -s tests -v
```

## Output

Plasmid and feature modes write a CSV with one row per annotated feature. The output includes feature coordinates and length, status, covered bases, breadth, mean depth, mean identity, and an `evidence_reads` count for reads that meet the per-read feature criteria. Coordinates in the result are 1-based and inclusive. `EVIDENCE_DETECTED` requires the feature to pass the configured breadth, mean-depth, and read-count thresholds; reads are also filtered using the alignment score and identity settings.

Evidence for separate features does not by itself show that they occur together on the same plasmid molecule or establish their order and orientation in an unknown assembly. PlasmidProbe reports evidence relative to the supplied reference and annotations.

## License

PlasmidProbe is released under the MIT License. Replace the copyright-holder placeholder in `LICENSE` before publishing.
