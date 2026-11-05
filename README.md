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
development build**, revision `553de47` (nightly built November 4, 2026).
The current wrapper requires this newer CLI: selecting `842adef` via
`--openms_container` alone fails because four MS1 options are unsupported there.
Historical comparisons use the retained
[old-build configuration overlay](results/end_to_end_comparison_20261104/old.config).

Build the two helper images from the repository root:

```bash
docker build -f containers/sdrf-pipelines.Dockerfile \
    -t openms-silac/sdrf-pipelines:0.1.7-ontology1 containers

docker build -f containers/python.Dockerfile \
    -t openms-silac/python:3.11-procps containers
```

The SDRF recipe uses the digest-pinned public SDRF-pipelines 0.1.7 image and
downloads its ontology cache during the build. The parser uses that cache for
offline ontology validation. Building needs network access; runtime parsing
does not need it for the validated fixtures. The cache is readable by non-root
users under `/opt/sdrf-cache`.

Select these helpers with `--sdrf_container` and `--python_container`; their
defaults are the build tags above. They are local images and must be transferred
or published before use on another machine. Enable containers with `-profile docker`.

QPX enrichment and pMultiQC reporting are enabled by default. They require the
**locally patched consumer images** selected by `--qpx_container` and
`--pmultiqc_container` in [nextflow.config](nextflow.config). The two helper
builds above do not provide them. The default `sha256:…` references identify
images on the tested Docker daemon, not published registry images. Follow the local
[consumer build instructions](containers/CONSUMERS.md) and
[build manifest](containers/consumer-builds.json), or transfer the verified image
with Docker save/load. Verify the resulting image ID; a rebuild may require
selecting the new verified ID through the corresponding container parameter.
Consumer recipes and patches under `containers/` are Git-ignored and must be
transferred separately when setting up another checkout. An unpatched released
QPX image does not provide the tested SILAC enrichment behavior.

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

Outputs are published directly under `--outdir`. Use a separate output
directory for each study, for example `--outdir results/my_study`.

```text
OUTDIR/
├── sdrf/
│   ├── study.sdrf.tsv
│   ├── study.sdrf_openms_design.tsv
│   ├── study.sdrf_config.tsv
│   ├── silac_chemistry.json
│   └── silac_config_<run-id>.json
├── quant_tables/
│   ├── experiment.mzTab
│   └── experiment.consensusXML
├── qpx/                         # enabled by default
└── pmultiqc/                    # enabled by default
    ├── *.html
    └── *_data/
```

With `--enable_qpx true`, the workflow enriches native feature, PSM, and
protein-group Parquet tables with SDRF metadata and publishes the QPX dataset.
Strict QPX schema validation is a separate check. Keep the mzTab and
consensusXML companions for peptide-level targeted/requantified evidence that
the QPX feature table does not preserve.
The report consumes enriched QPX and matching design/SDRF/mzTab context.
[assets/multiqc_config.yml](assets/multiqc_config.yml) supplies report identity
and file-recognition patterns; no quantms logo is staged. Its title also
determines the HTML/data basename in the tested MultiQC version. Report plots
summarize the available evidence; they do not establish native-ratio accuracy.

Use `--enable_pmultiqc false` for QPX without reporting. To disable both stages,
pass `--enable_qpx false --enable_pmultiqc false`. Reporting requires QPX.
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
| `--silac_join_channels` | `true` |
| `--silac_targeted_quantification` | `false` |
| `--silac_requantify_min_correlation` | `0.5` |
| `--silac_requantify_correlation_ratio` | `0.3` |
| `--silac_id_mz_reference` | `peptide` |
| `--silac_match_between_runs` | `false` |
| `--openms_peakpicking` | `false` |
| `--enable_qpx` / `--enable_pmultiqc` | `true` / `true` |

Targeted quantification is disabled by default. Enable it explicitly with
`--silac_targeted_quantification true` to extract identified peptides without
a detected multiplet. Retain mzTab and consensusXML: they preserve peptide-level
targeted and re-quantified evidence counts that QPX feature tables do not carry.
QPX protein-group tables retain the corresponding protein-level counts.

The `low_memory` profile runs one task at a time; it does not cap per-task RAM.
The `test_duplex` and `test_triplex` profiles require local fixtures and data
that are absent from a fresh clone.

## Local analysis and maintenance

Published mzTab and consensusXML retain native quantities and evidence for
downstream analysis. Local analysis prototypes are retained under
`bin/legacy/`; compatibility with the current output layout remains unverified.

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

Historical multithreaded runs on earlier OpenMS builds showed small peptide-ratio
differences. The nightly dependency evaluation used retained IDs and one thread;
it does not establish end-to-end identification equivalence for the new build.
A completed old/new comparison from fresh searches covered E. coli duplex and
triplex, with all 56 tasks successful. Identification populations differ slightly
after rescoring. A separate duplex check established bounded one-/two-thread
replay; this does not establish reproducibility across all datasets.
Shared Percolator input preparation changed despite unchanged
external search/rescoring binaries. OpenMS also records a hardcoded Percolator
version of 3.07; the pinned image contains Percolator 3.06.0.
