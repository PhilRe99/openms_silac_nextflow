# Comet fragment-tolerance resolution

`CREATE_INPUT_CHANNEL` reads per-run metadata. `RESOLVE_SILAC_CONFIG` then
resolves and validates the Comet fragment bin tolerance before any spectrum
processing or decoy generation starts. Comet consumes the resolved value;
it does not choose between metadata and parameters again.

## Precedence and validation

1. If `--comet_fragment_mass_tolerance` is supplied, use that value in **Da**.
2. Otherwise, use `FragmentMassTolerance` and `FragmentMassToleranceUnit`
   from the converted SDRF/OpenMS config.
3. The effective unit must be `Da`, and the value must be finite, numeric,
   and at least `0.0001`. Invalid overrides fail instead of falling back.

Metadata in ppm requires an explicit Da override. There is no fixed
ppm-to-Da conversion independent of fragment m/z. Selecting a Comet bin
tolerance is an analysis choice, not a claim of equivalence to a metadata
tolerance in ppm.

Correct inaccurate metadata in the SDRF. For a deliberate analysis override,
add this option to the normal processing command:

```bash
--comet_fragment_mass_tolerance 0.02
```

The same existing parameter can be supplied through a Nextflow parameter file.
`explicit_override` means the parameter supplied the value; it does not identify
whether that parameter came from the CLI, a parameter file, or custom config.

## Published resolution output

Each successful resolver task publishes:

```text
<outdir>/<input SDRF name without .sdrf.tsv>/configuration/
  silac_config_<run-id>.json
```

The JSON retains the existing SILAC chemistry, variable/binary modifications,
original OpenMS row, and matching experimental-design rows. It adds:

```json
{
  "comet_fragment_mass_tolerance": {
    "value": 0.02,
    "unit": "Da",
    "source": "explicit_override",
    "metadata_value": "20",
    "metadata_unit": "ppm"
  }
}
```

With no override, the source is `metadata`. Original metadata fields remain
unchanged; resolved tolerance, unit, and source are attached separately to
the run metadata used by Comet. Numeric formatting may be canonicalized
(for example, `0.020` becomes `0.02`); the numeric tolerance is preserved.

These outputs describe **resolved SILAC configuration and Comet fragment
settings**, not every effective setting or implicit tool default. They can
appear before the complete experiment finishes and do not certify success.

To check execution, inspect the matching Comet task's `work/<hash>/.command.sh`.
Its `-fragment_mass_tolerance` and `-fragment_error_units` must match the JSON's
`value` and `unit`. Intermediate spectra, IDs, and command files still remain
in the Nextflow work directory.

## Helper compatibility and tests

The helper still accepts its five positional arguments. An optional
`--comet-fragment-mass-tolerance` argument supplies the override directly
(the helper uses hyphens; the workflow parameter uses underscores).

Run the helper and baseline tests from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

All runs must resolve successfully and pass the existing homogeneous-chemistry
check before downstream processing starts. Other search settings,
quantification settings, containers, and defaults are unchanged by this increment.
