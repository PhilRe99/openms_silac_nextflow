# OpenMS SILAC Nextflow prototype

This repository incrementally migrates an existing Bash workflow for
MS1-labeled SILAC processing with OpenMS to Nextflow DSL2.

The project is intentionally small. Each workflow step is added and checked
against the Bash reference before the next step is introduced. It is not yet a
complete SILAC analysis pipeline.

## Current scope

The current Phase 1 workflow:

1. validates an SDRF file;
2. converts it to `openms.tsv` and `experimental_design.tsv` with
   `sdrf-pipelines`;
3. creates one `[meta, mzML]` Nextflow record per physical input file; and
4. keeps the OpenMS configuration and experimental design available as
   separate outputs for later stages.

Identification, spectrum preparation, MS1 quantification, and label-chemistry
cohorting have not yet been implemented in Nextflow.

## Requirements

- Nextflow 25.10.0 or newer
- A local environment containing the `parse_sdrf` command
- An SDRF metadata file
- The corresponding mzML files in a local directory

The current development setup uses the `silac` Conda environment from the
companion Bash project.

## Run the Phase 1 example

With `openms_silac_bash` and this repository as sibling directories:

```bash
conda activate silac
cd ~/thesis/openms_silac_nextflow

nextflow run . \
    --input ../openms_silac_bash/metadata/PXD003327/PXD003327.sdrf.tsv \
    --mzml_dir ../openms_silac_bash/data/PXD003327/mzml
```

After a code change, reuse completed work with:

```bash
nextflow run . \
    --input ../openms_silac_bash/metadata/PXD003327/PXD003327.sdrf.tsv \
    --mzml_dir ../openms_silac_bash/data/PXD003327/mzml \
    -resume
```

The workflow should print three run records for PXD003327:

```text
[[id:Chris_Ecoli_4-1], <path>/Chris_Ecoli_4-1.mzML]
[[id:Chris_Ecoli_1-1], <path>/Chris_Ecoli_1-1.mzML]
[[id:Chris_Ecoli_1-2-4], <path>/Chris_Ecoli_1-2-4.mzML]
```

The record order is not part of the workflow contract. Verify that all three
paths exist and that the generated OpenMS tables contain three unique physical
runs without duplicates caused by SILAC channel rows.

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
