# OpenMS SILAC Nextflow prototype

This repository incrementally migrates an existing Bash workflow for
MS1-labeled SILAC processing with OpenMS to Nextflow DSL2.

The project is intentionally small. Each workflow step is added and checked
against the Bash reference before the next step is introduced. It is not yet a
complete SILAC analysis pipeline.

## Current scope

The workflow parses an SDRF, resolves the SILAC configuration for each mzML,
performs OpenMS identification, and runs `MS1LabeledWorkflow` once for the
experiment. Each invocation must use one compatible SILAC label chemistry, so
the duplex and triplex subsets of PXD003327 run separately.

## Requirements

- Nextflow 25.10.0 or newer
- Docker for the pinned `sdrf-pipelines` container
- The local OpenMS and Python tools used by the remaining processes (the
  development setup uses the `silac` Conda environment)
- The sibling `openms_silac_bash` repository with the PXD003327 mzML files and
  `db/ecoli_k12_uniprot_reviewed.fasta`

## Run the PXD003327 test profiles

From this repository:

```bash
conda activate silac
nextflow run . -profile docker,test_duplex -resume
nextflow run . -profile docker,test_triplex -resume
```

The `test_duplex` profile selects `PXD003327_duplex.sdrf.tsv` (two mzML runs),
and `test_triplex` selects `PXD003327_triplex.sdrf.tsv` (one mzML run). Both use
the PXD003327 mzML directory and E. coli target FASTA in the sibling Bash
repository. The `docker` profile enables the SDRF parser container; the test
profiles select input data. `-resume` reuses completed tasks when their inputs
and settings match.

If the companion repository is elsewhere, override either path on the command
line, for example:

```bash
nextflow run . -profile docker,test_duplex \
    --mzml_dir /path/to/PXD003327/mzml \
    --database /path/to/ecoli_k12_uniprot_reviewed.fasta \
    -resume
```

Completed results are published under `results/PXD003327_duplex/` and
`results/PXD003327_triplex/`. Inspect `experiment.mzTab`,
`experiment.consensusXML`, and the generated SDRF configuration and design
TSV files in each directory.

## Check and export SILAC ratios

First check an existing duplex result against the known PXD003327 peptide
ratios:

```bash
python3 bin/analyze_ms1_mztab.py \
    results/PXD003327_duplex/experiment.mzTab \
    results/PXD003327_duplex/PXD003327_duplex.sdrf_openms_design.tsv \
    --benchmark pxd003327
```

The preserved checkpoint in `tests/fixtures/pxd003327_duplex_baseline.json`
records the command, input SHA-256 checksums, source software version, row
counts, and benchmark medians from the completed duplex result. Its software
version is `3.6.0-pre-HEAD-2026-09-30`. This is a migration reference, not a
universal acceptance threshold for other datasets or OpenMS versions.

Export reusable ratio tables from the same scientific outputs:

```bash
python3 bin/export_silac_ratios.py \
    --mztab results/PXD003327_duplex/experiment.mzTab \
    --consensusxml results/PXD003327_duplex/experiment.consensusXML \
    --design results/PXD003327_duplex/PXD003327_duplex.sdrf_openms_design.tsv \
    --outdir results/PXD003327_duplex/tables
```

The exporter uses Python's standard library and writes `peptide_ratios.tsv`
and `protein_ratios.tsv`. Use matching mzTab, consensusXML, and design files
from one MS1LabeledWorkflow invocation. The executed label definitions,
reference channel, and normalization setting are read from consensusXML's
MS1LabeledWorkflow processing record, so the original resolver JSON files
are not required. Sample IDs and condition text come directly from the design;
the exporter does not infer biological conditions from isotope labels.

Each peptide row represents `(sequence, modifications, fraction_group,
numerator_channel, reference_channel, normalization)`. Repeated exports across
charges or protein accessions are deduplicated; conflicting stored ratios or
support counts cause an error. The `peptidoform` column retains OpenMS's
modified-sequence annotation; it is not converted to ProForma. All non-decoy
peptide rows are eligible, including peptides with multiple K/R residues.

Each protein row represents `(full protein-group membership, fraction_group,
numerator_channel, reference_channel, normalization)`. Group membership is
sorted and semicolon-separated. Protein quantities come from the stored
`fraction_group_level_ratio` fields, never from quotients of protein abundance
columns. Peptide quantities come from mzTab's stored `MS1Label:peptide_ratio`
fields, not per-feature evidence ratios.

Both tables include fraction-ordered run names, channel definitions, sample
IDs, conditions, ratio, log2 ratio, and OpenMS's stored `ratio_count`. Raw and
normalized estimates are separate rows when ratio normalization was enabled.
Reference-to-reference comparisons are omitted. Peptide support counts refer
to OpenMS's contributing evidence counts; protein support counts refer to its
contributing peptide-ratio counts, not independent biological replicates.

Every exported entity has a row for each design fraction group and each
non-reference channel. Missing estimates have `null` quantities and
`ratio_status=not_reported`; this does not distinguish undetected signal from
insufficient evidence. Non-positive ratios retain their stored value with no
log2 value; non-finite ratios use `null` and `ratio_status=non_finite`. Missing
values are never replaced by zero. The current duplex checkpoint has 954
reported peptide ratios across 701 modified peptides, and 182 reported protein
ratios across 117 groups. Total table rows also include missing comparisons.

Run the small format tests and, when the local duplex outputs are available,
the saved-baseline regression:

```bash
python3 -m unittest discover -s tests -v
```

These helpers run independently of Nextflow.

## Native QPX output

Add `--enable_qpx` to export QPX directly from `MS1LabeledWorkflow` using
the configured OpenMS container. For the duplex fixture moved into `testdata/`:

```bash
nextflow run . -profile docker,test_duplex \
    --input testdata/PXD003327_duplex.sdrf.tsv \
    --enable_qpx -resume
```

The MS1 task adds `-out_qpx qpx` and publishes the directory under
`results/PXD003327_duplex/qpx/`. OpenMS writes `quantms.feature.parquet`,
`quantms.psm.parquet`, and `quantms.pg.parquet` alongside the existing mzTab
and consensusXML outputs. The files are also emitted as
`MS1_LABELED_WORKFLOW.out.qpx`. The default `enable_qpx = false` omits this
optional output. There is no separate `qpxc` conversion task on this path.

Native export does not run the QPX toolkit's strict schema validator or SDRF
metadata enrichment. Keep the explicit SILAC ratio tables for checking stored
ratio estimates; enabling export alone does not establish compatibility with
a particular downstream QPX reader.

## Repository layout

```text
main.nf                                      Pipeline entry point
workflows/silac.nf                           High-level workflow orchestration
modules/local/sdrf_parsing/main.nf           SDRF validation and conversion
subworkflows/local/create_input_channel/     Per-file channel construction
bin/                                         Validated helpers for later phases
bash_bp/                                     Bash regression reference
```

See `AGENT.md` for the development principles, scientific context, success
criteria, and incremental roadmap.
