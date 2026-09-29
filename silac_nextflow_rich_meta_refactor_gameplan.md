# Gameplan: Refactor SILAC Nextflow Prototype to Rich `meta` Records

## Goal

Refactor the current Nextflow prototype so that **per-run processing metadata lives in the Nextflow `meta` map**, similar to quantms.

Do **not** redesign the workflow or change scientific behavior.

The current downstream pattern is roughly:

```text
[meta, mzML, silac_config]
        ↓
PeakPicker
        ↓
[meta, mzML, silac_config]
        ↓
Comet
        ↓
[meta, mzML, idXML, silac_config]
        ↓
PeptideIndexer
        ↓
PSMFeatureExtractor
        ↓
Percolator
```

Change it to:

```text
[rich_meta, mzML]
        ↓
PeakPicker
        ↓
[rich_meta, mzML]
        ↓
Comet
        ↓
[rich_meta, mzML, idXML]
        ↓
PeptideIndexer
        ↓
PSMFeatureExtractor
        ↓
Percolator
```

`experimental_design.tsv` must remain a separate experiment-level artifact/channel. Do **not** put the whole design into `meta`.

The existing Python SILAC resolver remains the source of truth for SILAC chemistry. Do **not** rewrite its biological logic in Groovy.

---

## 1. Preserve the Current Scientific Behavior

Do not change:

- SDRF parsing
- `sdrf-pipelines`
- SILAC label interpretation
- Comet parameters
- modification chemistry
- binary modification grouping
- target-decoy generation
- PeakPicker behavior
- PeptideIndexer parameters
- PSMFeatureExtractor parameters
- Percolator parameters
- FDR thresholds
- database handling
- current container choices unless required by an existing module

The Bash/prototype behavior remains the regression reference.

Only change **where metadata is carried**.

---

## 2. Enrich `meta` When Creating the Initial Run Channel

Currently `CREATE_INPUT_CHANNEL` likely creates something close to:

```nextflow
[
    [id: run_id],
    mzml
]
```

Expand that `meta` map with the useful generic fields already present in the SDRF-generated OpenMS config table.

Use names that are simple and consistent.

Target conceptually:

```nextflow
meta = [
    id:                           run_id,
    experiment_id:                experiment_id,

    acquisition_method:           ...,
    label_type:                   ...,
    enzyme:                       ...,

    fixed_modifications:          ...,
    base_variable_modifications:  ...,

    precursor_mass_tolerance:       ...,
    precursor_mass_tolerance_unit:  ...,
    fragment_mass_tolerance:        ...,
    fragment_mass_tolerance_unit:   ...,

    dissociation_method:          ...
]
```

Then emit:

```text
[meta, mzML]
```

Do not force every OpenMS TSV column into `meta`. Only include values already relevant to the current workflow/search.

The key architecture is:

```text
per-run processing metadata
        ↓
      meta

experiment-level structure
        ↓
experimental_design.tsv
```

---

## 3. Keep `RESOLVE_SILAC_CONFIG` Unchanged Scientifically

The existing resolver process should still call the validated Python helper and produce:

```text
silac_config_<run>.json
```

The Python helper remains responsible for deriving:

```text
labels
label_channels
label_modifications
variable_modifications
binary_modifications
ffm_labels
```

Do not reproduce these calculations in Nextflow.

The resolver remains the SILAC-specific metadata extension point.

---

## 4. Immediately Convert Resolver Output Into Rich `meta`

Immediately after `RESOLVE_SILAC_CONFIG`, read its JSON once and merge the SILAC-specific values into the existing `meta` map.

Do this with a short `map` operator.

Do **not** create a new process or new subworkflow just for this.

Conceptually:

```nextflow
import groovy.json.JsonSlurper

resolved_runs = RESOLVE_SILAC_CONFIG.out.runs.map {
    meta, mzml, config ->

    def cfg = new JsonSlurper().parse(config.toFile())

    def resolved_meta = meta + [
        labels:                  cfg.labels,
        label_channels:          cfg.label_channels,
        label_modifications:     cfg.label_modifications,

        variable_modifications:  cfg.variable_modifications,
        binary_modifications:    cfg.binary_modifications,

        silac_labels:            cfg.ffm_labels
    ]

    tuple(resolved_meta, mzml)
}
```

### Important naming decision

Use:

```text
silac_labels
```

as the clear Nextflow `meta` name for:

```text
cfg.ffm_labels
```

Examples:

```text
duplex:
[][Lys8,Arg10]

triplex:
[][Lys4,Arg6][Lys8,Arg10]
```

For PXD003327, the expected chemistry is:

```text
Chris_Ecoli_4-1
    → [][Lys8,Arg10]

Chris_Ecoli_1-1
    → [][Lys8,Arg10]

Chris_Ecoli_1-2-4
    → [][Lys4,Arg6][Lys8,Arg10]
```

Do **not** confuse this chemistry string with the global experimental-design label numbers `1/2/3`.

Those retain their biological meaning:

```text
1 = light
2 = medium
3 = heavy
```

---

## 5. Stop Carrying `silac_config.json` Downstream

Once the JSON has been merged into `meta`, downstream tuples should no longer contain `config`.

Change:

```text
[meta, mzML, config]
```

to:

```text
[meta, mzML]
```

and:

```text
[meta, mzML, idXML, config]
```

to:

```text
[meta, mzML, idXML]
```

The JSON can still remain as a resolver output/work artifact for debugging or provenance.

Do not delete the resolver or JSON generation.

Just stop using the JSON as workflow state after meta enrichment.

---

## 6. Update `PEAK_PICKER`

Change only its tuple interface.

Before:

```nextflow
input:
tuple val(meta), path(mzml), path(config)

output:
tuple val(meta), path(prepared_mzml), path(config)
```

After:

```nextflow
input:
tuple val(meta), path(mzml)

output:
tuple val(meta), path(prepared_mzml)
```

Do not change the `PeakPickerHiRes` command or its scientific parameters.

---

## 7. Update `COMET`

This is the main functional refactor.

Currently Comet probably receives:

```text
meta
mzML
silac_config.json
database
```

and extracts values using `jq`.

Remove the config-file input.

Comet should instead use fields from `meta`.

Its input should conceptually become:

```text
[meta, mzML]
+
target_decoy_database
```

or the equivalent combined tuple already used by the workflow.

Use these values from `meta`:

```text
meta.enzyme
meta.fixed_modifications
meta.variable_modifications
meta.binary_modifications

meta.precursor_mass_tolerance
meta.precursor_mass_tolerance_unit

meta.fragment_mass_tolerance
meta.fragment_mass_tolerance_unit
```

Keep all existing Comet flags exactly the same:

```text
missed cleavages
instrument
isotope_error
reindex
force
threads
etc.
```

### Lists

Ensure list-valued metadata becomes valid command-line argument sequences.

Conceptually:

```nextflow
def variable_mods = meta.variable_modifications
    .collect { "\"${it}\"" }
    .join(' ')

def binary_mods = meta.binary_modifications
    .join(' ')
```

Use the simplest form compatible with the existing module.

Do not alter the actual modification contents.

---

## 8. Update Comet Output Tuple

Change:

```text
[meta, mzML, comet.idXML, config]
```

to:

```text
[meta, mzML, comet.idXML]
```

The enriched `meta` now carries all run-specific state.

---

## 9. Update `PEPTIDE_INDEXER`

Only change tuple shape.

Input:

```text
[meta, mzML, comet.idXML]
+
target_decoy_database
```

Output:

```text
[meta, mzML, indexed.idXML]
```

Keep:

```text
-decoy_string DECOY_
-decoy_string_position prefix
```

unchanged.

---

## 10. Update `PSM_FEATURE_EXTRACTOR`

Input:

```text
[meta, mzML, indexed.idXML]
```

Output:

```text
[meta, mzML, features.idXML]
```

Do not change its command.

---

## 11. Update `PERCOLATOR`

Input:

```text
[meta, mzML, features.idXML]
```

Output:

```text
[meta, mzML, percolator.idXML]
```

Keep the current per-file Percolator behavior and parameters unchanged:

```text
-use_subprocess true
-score_type pep
-score:fdr ...
```

---

## 12. Do Not Implement Grouping Yet

Do **not** add:

- `groupTuple`
- chemistry cohorts
- `MS1LabeledWorkflow`
- group-specific experimental designs
- new manifests
- another metadata-preparation subworkflow

This task ends when Percolator produces:

```text
[meta, mzML, idXML]
```

where `meta` already contains:

```text
id
generic search metadata
SILAC search metadata
silac_labels
```

Later we should be able to group directly with something like:

```nextflow
PERCOLATOR.out.runs
    .map { meta, mzml, idxml ->
        tuple(meta.silac_labels, meta, mzml, idxml)
    }
    .groupTuple(by: 0)
```

But **do not implement that in this task**.

---

## 13. Desired Final Dataflow After This Refactor

```text
SDRF
 │
 ▼
SDRF_PARSING
 │
 ├── openms.tsv
 └── experimental_design.tsv
          │
          │
CREATE_INPUT_CHANNEL
          │
          ▼
[generic-rich-meta, mzML]
          │
          ▼
RESOLVE_SILAC_CONFIG
          │
          ├── silac_config.json
          │
          ▼
merge SILAC config into meta
          │
          ▼
[fully-resolved-meta, mzML]
          │
          ▼
PEAK_PICKER
          │
          ▼
[meta, prepared_mzML]
          │
          ▼
COMET
          │
          ▼
[meta, mzML, comet.idXML]
          │
          ▼
PEPTIDE_INDEXER
          │
          ▼
[meta, mzML, indexed.idXML]
          │
          ▼
PSM_FEATURE_EXTRACTOR
          │
          ▼
[meta, mzML, features.idXML]
          │
          ▼
PERCOLATOR
          │
          ▼
[meta, mzML, percolator.idXML]
```

---

## 14. Acceptance Criteria

After the refactor, run the existing PXD003327 test with `-resume`.

Verify:

1. The workflow still produces exactly three run records.
2. The same mzMLs are processed.
3. Comet receives the same enzyme, tolerances, fixed modifications, variable modifications and binary groups as before.
4. The generated Comet command is scientifically equivalent to the pre-refactor command.
5. PeptideIndexer, PSMFeatureExtractor and Percolator still succeed unchanged.
6. `silac_config.json` is no longer present in tuples after the enrichment map.
7. Final Percolator records have this shape:

```text
[meta, mzML, idXML]
```

8. Inspecting `meta` for PXD003327 shows approximately:

```text
Chris_Ecoli_4-1
    silac_labels = [][Lys8,Arg10]

Chris_Ecoli_1-1
    silac_labels = [][Lys8,Arg10]

Chris_Ecoli_1-2-4
    silac_labels = [][Lys4,Arg6][Lys8,Arg10]
```

9. `experimental_design.tsv` is still preserved separately.
10. No MS1 grouping or cohort logic has been added yet.

---

## Scope Guard for the Coding Agent

> Keep this refactor minimal. The purpose is only to move stable per-run processing state from `silac_config.json` into Nextflow `meta`, closer to the quantms architecture. Do not redesign unrelated modules, introduce new helper layers, change scientific parameters, add validation frameworks, or proceed to MS1 grouping. Reuse the existing resolver and existing working tool commands.
