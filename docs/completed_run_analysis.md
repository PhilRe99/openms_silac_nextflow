# Completed-run scientific analysis

Run from any directory using Python 3.10 or newer:

```bash
python3 /path/to/openms_silac_nextflow/bin/analyze_completed_run.py /path/to/results
```

The default output is `RESULTS/analysis/report.html`, with JSON metrics and a
separate folder of TSV tables per mzTab. Existing scientific inputs are read
without modification. A parent result directory can contain multiple runs and
cohorts; these are analyzed separately. Re-running replaces analysis outputs;
use a fresh `--outdir` when preserving an earlier report.

```bash
# Current duplex, entirely without plotting dependencies
python3 bin/analyze_completed_run.py results/PXD003327_duplex --no-plots

# Current HeLa run, once quantification outputs have been published
python3 bin/analyze_completed_run.py \
    results/pxd039580_raw_all6_20261006_210459 --dataset PXD039580

# Cohort-local experimental_design.tsv files are discovered automatically
python3 bin/analyze_completed_run.py results/cohort_validation_mixed

# Explicit files for unconventional naming or ambiguous designs
python3 bin/analyze_completed_run.py /data/experiment.mzTab \
    --consensusxml /data/experiment.consensusXML \
    --design /data/design.tsv --outdir /data/analysis
```

The matching design is sought beside the mzTab, then in its ancestors up to
the requested result directory. Its run/fraction/channel map must match the
consensusXML; two different matching designs require `--design`. Explicit
file overrides are limited to one mzTab. Keep each pipeline invocation in a
distinct published directory.

The script targets this repository's OpenMS outputs, rather than arbitrary
proteomics formats. Without consensusXML or a matching design it retains
available mzTab identification summaries and explains why quantitative
analysis is unavailable. With no mzTab it produces a diagnostic report and
summarizes available execution traces. Recorded trace tasks may omit work
still running or pending. A malformed bundle does not prevent analysis of
other bundles. Invalid CLI/configuration inputs exit 1; ordinary missing
analyses are reported with exit 0, or exit 2 using `--strict`.

Matplotlib produces PNG and PDF plots when available; otherwise a diagnostic
is added. `--no-plots` disables plotting. SciPy is optional and used only for
configured biological tests. The script does not install packages or access
the network at runtime.

## Dataset questions and evidence

Definitions live in `conf/analysis_profiles.json`. Exact run stems map to
mixture amounts and replicate identities. Dataset detection uses accession
text in the input path, or explicitly permitted known run sets. Use
`--dataset none` to disable expectations or `--dataset ACCESSION` to override
detection. Every profile and its source links are preserved in JSON metrics.

| Dataset | Scientific question | Interpretation |
| --- | --- | --- |
| PXD003327 | Does quantification recover the known isotope mixtures for single-K and single-R peptides? | Expectations are residue-specific; pooled protein ratios have no universal ground truth. |
| PXD039580 | How accurately and reproducibly does DDA quantify the two HeLa mixtures? | Mix1 L:M:H = 70:15:15; Mix2 = 20:40:40. Each has three technical injections. These are not biological treatment replicates. |
| MSV000079761 | Does the yeast mixture recover H/L = 1? | One mixture does not establish replicate precision. Accession detection or explicit selection is required because `1to1` is a generic filename. |

PXD003327 uses Table I of the original
[PyQuant publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC4974355/).
PXD039580 uses Fig. 4 and the sample-preparation methods of the
[QuaNPA publication](https://www.nature.com/articles/s41467-023-43919-3);
its run mapping follows the locally curated SDRF and author benchmark file
names. The benchmark cannot answer IFNγ treatment questions: that experiment
is deposited as PXD038915. A DDA-only analysis cannot establish a comparison
with DIA without matched external results. MSV000079761's 1:1 expectation is
documented in the original
[MassIVE deposit](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=6a13399e5a4145b2b87dd1b5a1d55948);
no associated peer-reviewed publication was established here.

## Tables and statistical definitions

| File in each bundle folder | Contents |
| --- | --- |
| `peptide_ratios.tsv`, `protein_ratios.tsv` | Validated stored OpenMS ratios with chemistry, conditions, support counts, and explicit missingness. |
| `summary.tsv` | Coverage, median ratios, inclusive log2 quartiles, and support, separated by fraction group and normalization. |
| `benchmark.tsv` | Expected ratio, median bias, median log2 error, median absolute log2 error, log2 IQR, and class coverage. |
| `replicate_precision.tsv` | Per-entity sample CV on linear ratios and SD on log2 ratios, including observed/planned replicate counts. |
| `replicate_correlations.tsv` | Pairwise log2 Pearson correlations and median absolute log2 differences among shared quantified entities. |
| `biological_effects.tsv` | Configured paired protein-group effects with biological sample counts and optional p/q-values. |
| `metrics.json` | Findings, diagnostics, executed OpenMS settings, software annotations, SHA-256 input provenance, and the scientific profile. |

Coverage is measured against experiment-wide exported entities, including
those without a ratio in a particular group. Missing estimates are never
imputed or treated as biological absence. Raw and normalized results remain
separate. Linear medians average the two central linear ratios for even-sized
samples; log2 medians average in log space. Use raw ratios as the primary accuracy measure for deliberately
unequal mixtures; normalization can remove the intended displacement.

All pairwise channel comparisons are available for benchmarks and replicate
analysis. The executed reference is read from consensusXML, including
non-light references. H/M calculated as `(H/L)/(M/L)` is marked
`derived_from_stored_ratios`; it requires both source ratios. Dividing
aggregated protein ratios can differ from native OpenMS H/M quantification
with a medium reference. It is an explicitly derived estimate. The stored
ratio exports remain in their original reference orientation.

Precision summaries distinguish complete-case CV from CV using any two or
more available replicates. Planned but absent injections remain in the
completeness denominator. Fractions stay within their experimental-design
fraction groups. Duplicate groups assigned to one replicate are rejected;
combine fractions in the design or correct replicate identities. Separate
pipeline invocations are never pooled automatically. A CV among technical
injections measures technical precision. Correlation requires at least three
shared entities and nonconstant vectors; it does not itself establish
accuracy or agreement.

## Configure another study or a biological contrast

Supply your own JSON with `--profiles my_study.json --dataset my_study`.
Amounts are optional; unknown studies can use replicate definitions or
condition contrasts without asserting quantitative ground truth. Keys in
`runs` must match spectral-file stems, excluding `.raw` or `.mzML`.
Channel keys are executed label definitions, rather than condition text or
assumptions about channel numbers. In a fractionated group, all member runs
must have identical scientific metadata.

This example defines two biological replicates and a label swap; add a third
independent replicate before a statistical test is available:

```json
{
  "schema_version": 1,
  "datasets": {
    "my_study": {
      "sources": ["https://example.org/study"],
      "peptide_classes": {},
      "runs": {
        "forward": {
          "biological_replicate": "1",
          "comparison": "drug_24h",
          "conditions": {"[]": "vehicle", "[Lys8,Arg10]": "drug"}
        },
        "reverse": {
          "biological_replicate": "2",
          "comparison": "drug_24h",
          "conditions": {"[]": "drug", "[Lys8,Arg10]": "vehicle"}
        }
      },
      "contrasts": [
        {"name": "drug_vs_vehicle_24h", "numerator": "drug",
         "denominator": "vehicle", "comparison": "drug_24h"}
      ]
    }
  }
}
```

Repeated injections receive the same `biological_replicate`; their log2
effects collapse by median before testing. Do not assign independent
biological IDs to repeated injections of one sample. To request technical
precision as well, add `replicate_set`, `replicate_type: "technical"`, and a
distinct `replicate` identifier to each injection. Keep dose, time, batch,
and distinct treatment questions in separate `comparison` values.

For a known homogeneous mixture, define `"peptide_classes": {"all": {}}`
and add e.g. `"amounts": {"all": {"[]": 1, "[Lys8,Arg10]": 2}}` to each
applicable run. An empty residue selector means all peptides and proteins;
nonempty selectors use exact residue counts and apply only to peptides.

Biological contrasts compare two conditions present within the same mixture;
this does not fit arbitrary between-run designs or shared-reference models.
With at least three biological replicates and nonzero variance, optional
[SciPy one-sample t-tests](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_1samp.html)
test whether the mean log2 effect differs from zero. BH adjustment is applied
within each contrast and normalization in each bundle. Missing samples remain
missing. Zero-variance and under-replicated effects are descriptive only.
This assumes independent biological replicates and approximately normal
log2 effects; it is an exploratory analysis, without the moderated models,
covariates, pathway enrichment, or localization needed for some studies.
Protein groups retain full membership and must not be interpreted as a
unique gene when ambiguous. No universal significance or scientific
acceptance threshold is imposed.

## Verification

```bash
python3 -m unittest discover -s tests -v
```

Small fixtures cover channel reference changes, missing derived ratios,
published mixture expectations, planned-but-missing injections, label swaps,
technical-repeat collapse, BH adjustment, ambiguous designs, HTML escaping,
and incomplete outputs. Existing local duplex outputs additionally exercise
the preserved quantitative baseline. Real mixed-cohort outputs can be used
to check discovery and optional PNG/PDF generation.
