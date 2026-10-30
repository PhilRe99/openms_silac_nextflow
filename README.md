# OpenMS SILAC Nextflow prototype

Nextflow DSL2 workflow for DDA peptide identification and MS1 SILAC
quantification with OpenMS. It accepts SDRF metadata, Thermo RAW or mzML files,
and a target-only FASTA. Integration into quantms is a longer-term goal.

The workflow validates metadata, resolves SILAC chemistry, prepares spectra,
searches with Comet, indexes peptides, rescores with Percolator, and runs
`MS1LabeledWorkflow`. Identification uses native `.idparquet` bundles;
quantification produces mzTab and consensusXML, with optional QPX output.

## Setup

Requirements: Nextflow **25.10.0+**, a compatible Java runtime, and Docker.
OpenMS is pinned by digest in [nextflow.config](nextflow.config) to a **3.7
development build**, revision `842adef`.

Build the two helper images from the repository root:

```bash
docker build -f containers/sdrf-pipelines.Dockerfile \
    -t openms-silac/sdrf-pipelines:dea8ba9-procps containers

docker build -f containers/python.Dockerfile \
    -t openms-silac/python:3.11-procps containers
```

The SDRF recipe requires the existing base image
`openms-silac/sdrf-pipelines:dea8ba9`; this repository does not build that base.
Enable containers with `-profile docker`.

## Run

Supply an SDRF whose spectral filenames match the input directory and a
**target-only FASTA** for the organism. The workflow generates decoys.

For mzML inputs:

```bash
nextflow run . -profile docker \
    --input /path/to/study.sdrf.tsv \
    --mzml_dir /path/to/mzml \
    --database /path/to/targets.fasta \
    --outdir results/study \
    -resume
```

For Thermo RAW inputs:

```bash
nextflow run . -profile docker,low_memory \
    --input /path/to/study.sdrf.tsv \
    --spectra_dir /path/to/raw \
    --convert_raw \
    --database /path/to/targets.fasta \
    --outdir results/study_raw \
    -resume
```

`--spectra_dir` takes precedence over `--mzml_dir`. RAW conversion uses native
Thermo centroiding. For high-resolution profile mzML, add
`--openms_peakpicking true`; optionally restrict it with `--peakpicking_ms_levels 2`.
Already centroided mzML needs no additional picking. mzXML is unsupported.

Each invocation requires **one compatible SILAC chemistry**; process duplex
and triplex subsets separately. Chemistry definitions are in
[assets/silac_chemistry.json](assets/silac_chemistry.json).
If SDRF fragment tolerance is in ppm, choose and document an explicit
`--comet_fragment_mass_tolerance VALUE` in **Da**; no automatic conversion is used.

## Outputs and settings

Outputs are published under `OUTDIR/<SDRF name without .sdrf.tsv>/`:

```text
study.sdrf_config.tsv
study.sdrf_openms_design.tsv
configuration/silac_config_<run-id>.json
experiment.mzTab
experiment.consensusXML
qpx/                            # optional: --enable_qpx
```

QPX contains native feature, PSM, and protein-group Parquet tables. Export does
not automatically validate the QPX schema or enrich it with SDRF metadata.
Intermediate mzML and idparquet bundles stay in `work/`. Keep `work/` and
`.nextflow/` to reuse completed tasks with `-resume`.

Key defaults are in [nextflow.config](nextflow.config); OpenMS command options
are in [conf/modules.config](conf/modules.config).

| Parameter | Default |
| --- | --- |
| `--cpus` | `2` for the main search/quantification processes |
| `--psm_fdr` / `--protein_level_fdr_cutoff` | `0.01` / `0.01` |
| `--silac_reference_channel` | `1` (channels numbered from 1) |
| `--silac_min_ratio_count` | `2` |
| `--silac_normalize_ratios` | `false` |
| `--silac_knock_out` / `--silac_requantify` | `true` / `true` |
| `--silac_id_mz_reference` | `peptide` |
| `--silac_match_between_runs` | `false` |
| `--openms_peakpicking` / `--enable_qpx` | `false` / `false` |

The `low_memory` profile runs one task at a time; it does not cap per-task RAM.
The `test_duplex` and `test_triplex` profiles require local fixtures and data
that are absent from a fresh clone.

## Local analysis and maintenance

Where the local analysis helpers are available, generate a report from a
completed result directory:

```bash
python bin/analyze_completed_run.py results/study --dataset none --no-plots
```

This requires Python **3.10+** and the Git-ignored helper modules in
`bin/legacy/`. The report is written to `results/study/analysis/report.html`.
Omit `--dataset none` to use configured expectations from
[conf/analysis_profiles.json](conf/analysis_profiles.json). Matplotlib enables
plots; SciPy enables configured biological tests. Analysis reads native
mzTab/consensusXML/design files and preserves missing ratios.

Data, results, investigation reports, tests, and standalone analysis prototypes
remain local and Git-ignored. A fresh clone contains the processing workflow
and curated metadata, but lacks those local materials.

For failed tasks, inspect `.command.sh`, `.command.out`, and `.command.err` in
the task's `work/` directory. Preview cache cleanup with
`bash bin/cleanup_run.sh`; after stopping active runs, use `--delete` to remove
execution caches. Published results and inputs are retained; old runs can no
longer be resumed.

## Scientific scope

This is a prototype with dataset-specific validation. Missing ratios and
quantitative biases remain unresolved; reported coverage does not establish
accuracy. Keep the pinned scientific baseline; further attrition experiments
are paused by default. Phosphosite localization and quantms integration remain
separate work.
