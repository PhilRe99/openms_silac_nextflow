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
