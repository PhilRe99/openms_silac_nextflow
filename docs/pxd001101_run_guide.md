# PXD001101 investigation and pilot run guide

Reviewed 2026-10-07 against the local files, workflow source, pinned parser
container, deposited results, and literature.

**The corrupted sister-repository pilot was removed on 2026-10-07. With the corrected
pilot supplied here and an explicit Comet fragment setting, the existing
workflow can be used for an exploratory peptide-ratio run. Full RAW processing
and phosphopeptide quantification have not been tested in this review. Site
localization and site-level inference require additional work.**

## Storage and local inventory

Use this layout for new datasets:

```text
/mnt/d/thesis/data/<accession>/
  raw/          original vendor inputs
  reference/    target FASTA and database provenance
  mzml/         optional retained conversions
```

The annotation repository supplies curated metadata and supporting evidence;
keep large spectra out of it. Nextflow work and published results can remain
on the Linux filesystem. The current pipeline stages inputs and Docker mounts
them; no hard-coded RAW-storage change is needed.

| Item | Inspected state |
| --- | --- |
| `../nf_silac_data` | Contains PXD039580: six DDA RAWs (~8.8 GiB), curated SDRF, PRIDE inventory, reference MaxQuant XML/mzTab, FASTAs and checksum material. It has no Git repository at its root. |
| `../sdrf_annotations/datasets/PXD001101` | Publication/supplement, deposited MaxQuant ZIPs, complete RAW mapping, evidence notes and full curated SDRF. The corrupted pilot was removed. |
| `/mnt/d/thesis/data/PXD001101/raw` | Exactly `H1975_pST_BR1_F1.raw` through `F12.raw`: 12 files, 5,717,737,420 bytes (~5.33 GiB). |
| `/mnt/d/thesis/data/PXD001101/reference/human_uniprot_reviewed.fasta` | 20,416 entries; all inspected headers are human Swiss-Prot entries, no `DECOY`/`REV_` headers. 13,666,795 bytes. |
| `/mnt/d/thesis/data/PXD001101/PXD001101.sdrf.tsv` | Older 330-row annotation: structured labels, unknown technical index and older treatment/dissociation representation. Use the curated source instead. |
| `expected_raw_files.txt` on D: | Lists the full 110-file study, not just the downloaded pilot. |

All 12 RAW sizes and SHA-256 hashes of their first 65,536 bytes agree with the
annotation repository's remotely retrieved header evidence. This checks size
and header identity, not every byte of each download. The inspected PRIDE
inventory has empty checksum fields for the pilot files. Record whole-file
hashes locally when freezing the run; they establish future identity rather
than agreement with an unavailable depositor checksum.

The PXD039580 RAWs still exist under `nf_silac_data`; this review did not move
or delete them. Apply the D: convention to future downloads and migrate older
copies explicitly when convenient.

## Experimental design and biological question

PXD001101 is the dataset for [Zhang et al., Proteomics 2015,
doi:10.1002/pmic.201400315](https://pubmed.ncbi.nlm.nih.gov/25404012/): early
EGFR signaling responses to EGF and erlotinib in lung adenocarcinoma. H3255
is erlotinib-sensitive; H1975 carries L858R/T790M and is resistant. The
experiment uses phosphopeptide enrichment and triplex SILAC.

The detailed treatment, instrument, fraction and search reconstruction is
supported by the local [evidence notes](../../sdrf_annotations/datasets/PXD001101/evidence/experiment_design.md),
[source index](../../sdrf_annotations/datasets/PXD001101/evidence/sources.md)
and [unresolved issues](../../sdrf_annotations/datasets/PXD001101/notes/unresolved.md).
Those files preserve the publication and submitted design evidence, including
protocol contradictions. Their early research-phase stop statements precede
the later SDRF outputs.

| Channel | Chemistry | Condition | Useful comparison |
| --- | --- | --- | --- |
| L / 1 | Unlabeled K/R | Serum-starved control | Denominator for EGF response |
| M / 2 | Lys4 / Arg6 | EGF, 100 ng/mL, 3 min | M/L = EGF/control |
| H / 3 | Lys8 / Arg10 | Erlotinib, 100 nM, then EGF | H/M = erlotinib+EGF/EGF |

Erlotinib pretreatment is reported as both 60 and 120 minutes in the sources;
the annotation correctly leaves it unresolved. There is no erlotinib-only
channel and no reported label swap. Equal protein input was mixed 1:1:1;
individual phosphopeptide ratios have biological changes and are not a known
1:1 accuracy benchmark.

The curated full design contains 330 channel rows for 110 RAWs and eight
biological units: five H3255 experiments (73 files) and three H1975 experiments
(37 files). Preserve submitted fractions, including H1975 BR2 F13 and the
older H3255 fraction gaps. The downloaded pilot is one H1975 mixture with
12 fractions and three channels. It has **one biological replicate**; fractions
and isotope channels cannot supply replication for a biological significance
test. H1975 alone cannot reproduce the sensitive-versus-resistant comparison.

The pilot uses Orbitrap Elite HCD acquisition, high-resolution MS1/MS2,
basic reversed-phase fractionation and TiO2 phosphopeptide enrichment. The
older H3255 subset instead uses Orbitrap XL ion-trap MS/MS with multistage
activation and SCX. This analyzer split matters for subsequent full-study
searches.

## Readiness findings and required changes

1. **Regenerate the pilot without losing repeated headers.** The former
   `../sdrf_annotations/datasets/PXD001101/output/PXD001101_H1975_BR1_pilot.sdrf.tsv`
   contained eight identical acetylation entries per row and was deleted on
   2026-10-07 at the user's request. Its generated config
   has no fixed carbamidomethylation and no phosphorylation, oxidation or
   deamidation. A dictionary-based TSV round trip can produce this corruption
   because the eight physical modification columns have the same name;
   the generating command was not available to establish the exact cause.
   The full curated SDRF is correct. This repository supplies a
   [corrected pilot](../assets/metadata/PXD001101_H1975_BR1_pilot.sdrf.tsv)
   selected from that source with `csv.reader`/`csv.writer`, preserving every
   physical cell and all eight modification columns. Neither sister source
   nor the copy on D: was overwritten.
2. **Use the pinned pipeline parser.** Read-only inspection and actual
   conversion with `openms-silac/sdrf-pipelines:dea8ba9-procps` gave the full
   study 110 entries for each label and eight fraction groups. The corrected
   pilot gives 12 config rows, 36 design rows, labels 1/2/3 twelve times each,
   and one fraction group. This container includes the channel-handling
   fixes absent from the annotation review's earlier parser. The older
   `output/openms/` artifacts are not the pipeline's current conversion.
3. **Specify Comet fragment bins in Da.** The Elite SDRF says 20 ppm, which
   the module intentionally rejects. Start with
   `--comet_fragment_mass_tolerance 0.02`, the existing local high-resolution
   search choice. This is a Comet
   binning choice, not a conversion of 20 ppm. The adapter documents it as
   half the bin size. Assess search yield and mass-error behavior after the
   pilot before fixing it for comparative work.
4. **Choose the ratio reference deliberately.** The command below uses
   `--silac_reference_channel 2` for directly stored H/M drug-response ratios
   and L/M. Invert L/M for M/L EGF response. Using reference 1 instead stores
   M/L and H/L; the analysis tool can derive H/M, but division of aggregated
   protein ratios can differ from native H/M quantification.

The corrected parser output contains fixed `Carbamidomethyl (C)` and variable
`Acetyl (Protein N-term)`, `Deamidated (N/Q)`, `Oxidation (M)` and
`Phospho (S/T/Y)`. Resolution succeeds for all 12 files, with chemistry
`[][Lys4,Arg6][Lys8,Arg10]`, eleven variable-modification entries including
four label modifications, and binary groups `0,0,0,0,0,0,0,1,1,2,2`.

No pipeline logic change is required just to supply the D: paths or this
triplex chemistry. Two missed cleavages, 7 ppm precursor tolerance, high-res
Comet mode and existing FDR settings are appropriate starting settings to
retain. The installed CometAdapter defaults to at most five variable
modifications per peptide; isotope modifications also consume the modification
budget. Review this for highly phosphorylated/missed-cleavage peptides as a
documented sensitivity comparison. It was not changed here.

## Run the pilot

From the pipeline repository, using the corrected input supplied here:

```bash
cd /home/philipp/thesis/openms_silac_nextflow

run_dir=results/pxd001101_h1975_br1_pilot
mkdir -p "$run_dir"

nextflow -log "$run_dir/nextflow.log" run . \
    -profile docker,low_memory \
    --input assets/metadata/PXD001101_H1975_BR1_pilot.sdrf.tsv \
    --spectra_dir /mnt/d/thesis/data/PXD001101/raw \
    --convert_raw \
    --database /mnt/d/thesis/data/PXD001101/reference/human_uniprot_reviewed.fasta \
    --comet_fragment_mass_tolerance 0.02 \
    --silac_reference_channel 2 \
    --cpus 2 \
    --outdir "$run_dir" \
    -with-report "$run_dir/execution_report.html" \
    -with-trace "$run_dir/trace.tsv" \
    -with-timeline "$run_dir/timeline.html" \
    -resume
```

The local Nextflow version is 26.04.6. All four selected container images are
already present. The inspected machine has ~15 GiB RAM, ~822 GiB free on
Linux and ~1.4 TiB free on D:. `low_memory` serializes tasks; it does not cap
or lower the memory requirement of an individual Comet/MS1 task. The installed
CometAdapter processes up to 20,000 spectra per batch by default. Observe the
execution trace before increasing concurrency. No reliable runtime or peak-RAM
estimate exists for this pilot yet.

Keep the default Linux `work/` for this first run and retain it with
`.nextflow/` for resume. Converted spectra and identifications live there;
the pipeline does not automatically publish the conversions to D:/mzml.
Do not pass the full 110-file SDRF while only twelve RAWs are present: missing
files cause input construction to fail.

The expected published folder is:

```text
results/pxd001101_h1975_br1_pilot/PXD001101_H1975_BR1_pilot/
  PXD001101_H1975_BR1_pilot.sdrf_config.tsv
  PXD001101_H1975_BR1_pilot.sdrf_openms_design.tsv
  experiment.mzTab
  experiment.consensusXML
```

After successful quantification, generate descriptive QC:

```bash
python bin/analyze_completed_run.py \
    results/pxd001101_h1975_br1_pilot/PXD001101_H1975_BR1_pilot \
    --dataset none
```

Use `--no-plots` if Matplotlib is unavailable. There is currently no
PXD001101 entry in `conf/analysis_profiles.json`; `--dataset PXD001101` would
fail. Do not add known mixture expectations of one or interpret the generic
protein effect tables as phosphosite tests. Adding a biological profile alone
would not supply site localization or independent replicates.

Confirm successful completion of all twelve file-level paths and the final
MS1 task, correct phospho/search modifications in `.command.sh` and resolved
JSON, one fraction group with all fractions and three channels, reference
channel 2, and quantified modified-peptide coverage/missingness. A successful
execution report alone does not validate the scientific ratios.

## What further analysis needs

The human reviewed FASTA is structurally suitable as a target-only pilot
database. Its SHA-256 is
`dffdff3807ec6eb7e752035b5475aaefa26d1b36e513e3344198d521a2169c42`.
Its release/download provenance is not recorded next to the file. The original
study used RefSeq 49 plus contaminants, so this search will not reproduce
identification counts or protein groups exactly. Freeze the FASTA provenance
and use a documented contaminant set for a finalized reanalysis. Mapping the
deposited RefSeq identifiers to UniProt requires sequence/identifier care,
especially for isoforms and site coordinates.

Current processing has no phosphosite-localization step, site-level error
control, or validated site aggregation. PSM/protein FDR does not establish
correct phosphosite assignment. OpenMS [PhosphoScoring](https://openms.de/documentation/TOPP_PhosphoScoring.html)
offers Ascore localization, but adding it requires validation of the MS2,
identification, scoring and downstream quantification contracts. An Ascore is
not automatically the MaxQuant localization probability used in the paper.
Different phosphorylated peptides from one protein can change differently;
protein-level rollups of enriched peptides are not total-protein abundance or
site occupancy measurements.

For a full-study run, the workflow also needs analyzer-aware search choices:
73 Elite files carry 20 ppm FTMS annotation and 37 XL files carry 0.5 Da ITMS
annotation. The current global fragment override and global `--instrument`
cannot express appropriate different settings simultaneously. Separate
analyzer-specific invocations are possible, but separate inference/output
bundles complicate comparisons. For joint processing, add per-run Comet
settings with an explicit conversion/search policy and validate both branches.
The global SILAC chemistry is compatible across the study; this is an
acquisition-setting issue rather than a duplex/triplex issue.

The converter's sample table numbers `MSstats_BioReplicate` by source identity:
the pilot displays 1, 2 and 3 for the three conditions. These numbers do not
make three independent biological experiments. Preserve the SDRF biological
unit/fraction-group mapping explicitly in any later paired statistical design.

## Literature and available comparison results

- [Zhang et al. 2015, the original dataset paper](https://pubmed.ncbi.nlm.nih.gov/25404012/)
  is the primary biological reference. Its class-I sites use localization
  probability >0.75; the study reports weaker erlotinib suppression in H1975
  than H3255. Replicating that result needs localized sites and both cell lines.
- [Zhang et al. 2017, quantitative tyrosine phosphoproteomics](https://pubmed.ncbi.nlm.nih.gov/28331001/)
  and its [PXD004373 deposit](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD004373)
  are a closely related extension involving erlotinib/afatinib and sensitive/
  resistant lung adenocarcinoma models. It is a useful later comparison;
  tyrosine-specific enrichment differs from the TiO2 global phosphopeptide
  preparation here, so design/enrichment compatibility needs separate review.
- [Frankenfield et al. 2025, Benchmarking SILAC Proteomics Workflows and Data
  Analysis Platforms](https://pubmed.ncbi.nlm.nih.gov/40315959/)
  compares five packages across DDA/DIA and evaluates accuracy, precision,
  completeness and false discoveries. It supports retaining independent
  mixture benchmarks and checking against an additional analysis platform.
  It does not independently validate this prototype or phosphosite pipeline.

The deposited `H1975_MaxQuantresults.zip` already contains `proteinGroups.txt`,
`peptides.txt` and `PhosphoSTYSites.txt`, alongside design and parameter tables.
The H1975 site table has 6,559 rows and replicate-specific localization and
ratio columns. This is an unfiltered table row count, not the number of
validated sites. The H3255 ZIP and Supplementary Table 1 provide corresponding
reference comparisons. Neither ZIP contains `mqpar.xml`, `evidence.txt` or
`msms.txt`. Compare BR1-specific ratios when evaluating this pilot, rather
than ratios pooled across all three H1975 biological replicates.

## Verification and provenance

The corrected pilot is an exact physical-row subset of the full source with
SHA-256 `0dae6fa4249a6c363d5b048248630d2be2f9dc8528cf86a70b73cc3d2cb745cf`.
It was selected using `characteristics[biological unit] == H1975_pST_BR1`.
Every selected row is unchanged. It has 36 rows, twelve unique RAWs, three
source identities and eight distinct physical modification entries per row.
The corrected pilot's SHA-256 is
`f7dfa2bdc0102ba63124d7e326d95912a2db61ad0a97bf53413b6a9c5745198d`.

The pinned container's `validate-sdrf --template ms-proteomics` completed with
“Everything seems to be fine. Well done.” It also explicitly warned that
ontology dependencies are absent, so this is **structural validation only**.
The sister review separately records earlier ontology-enabled validation of
the full source. Conversion warned about the protein-N-terminal acetylation's
missing residue target and inferred its N-terminal position; generated search
modifications were inspected directly and correct.

Full/pilot conversions and all twelve chemistry-resolution JSONs are retained
locally under `results/pxd001101_readiness/`. The folder named
`pinned_parser_existing_pilot/` contains the corrected-pilot conversion after
an initial container-user permission failure; its name does not denote a
successful conversion of the corrupted sister pilot. No mass-spectrometry
workflow or new biological analysis was executed in this review.
