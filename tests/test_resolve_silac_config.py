"""Per-run design lookup and Comet fragment-tolerance resolution."""

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from resolve_silac_config_experiment import (
    read_experimental_design,
    resolve_comet_fragment_tolerance,
)


class FragmentToleranceTests(unittest.TestCase):
    def test_metadata_da_without_override(self):
        result = resolve_comet_fragment_tolerance("0.8", "Da")
        self.assertEqual(result, {
            "value": 0.8, "unit": "Da", "source": "metadata",
            "metadata_value": "0.8", "metadata_unit": "Da",
        })

    def test_override_takes_precedence_and_preserves_metadata(self):
        for metadata_value, metadata_unit in (("0.8", "Da"), ("20", "ppm"), (None, None)):
            with self.subTest(unit=metadata_unit):
                result = resolve_comet_fragment_tolerance(metadata_value, metadata_unit, "0.02")
                self.assertEqual(result["value"], 0.02)
                self.assertEqual(result["unit"], "Da")
                self.assertEqual(result["source"], "explicit_override")
                self.assertEqual(result["metadata_value"], metadata_value)
                self.assertEqual(result["metadata_unit"], metadata_unit)

    def test_ppm_or_missing_unit_without_override_fails(self):
        for unit in ("ppm", None, ""):
            with self.subTest(unit=unit):
                with self.assertRaisesRegex(ValueError, "--comet_fragment_mass_tolerance"):
                    resolve_comet_fragment_tolerance("20", unit)

    def test_missing_invalid_or_nonfinite_values_fail(self):
        for value in (None, "", "bad", "nan", "inf", "-inf", "1e999", -1, 0.00009, True, False):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "finite number >= 0.0001 Da"):
                    resolve_comet_fragment_tolerance(value, "Da")

    def test_invalid_override_never_falls_back_to_valid_metadata(self):
        for value in (0, "0", "", "bad", "nan", "inf", -1, False):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "finite number >= 0.0001 Da"):
                    resolve_comet_fragment_tolerance("0.8", "Da", value)

    def test_minimum_accepted_value(self):
        self.assertEqual(resolve_comet_fragment_tolerance("0.0001", "Da")["value"], 0.0001)


class ResolverCommandTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.openms = self.root / "openms.tsv"
        self.design = self.root / "design.tsv"
        self.chemistry = self.root / "chemistry.json"
        self.output = self.root / "silac_config_run.json"
        self.design.write_text(
            "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n"
            "1\t1\trun.mzML\t1\t1\n1\t1\trun.mzML\t2\t2\n\n"
            "Sample\tMSstats_Condition\n1\tlight\n2\theavy\n"
        )
        self.chemistry.write_text(json.dumps({"light": [], "heavy": ["Lys8", "Arg10"]}))

    def run_resolver(self, value, unit, override=None):
        with self.openms.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=[
                "Filename", "Label", "VariableModifications", "FragmentMassTolerance", "FragmentMassToleranceUnit",
            ])
            writer.writeheader()
            writer.writerow({"Filename": "run.raw", "Label": "SILAC", "VariableModifications": "Oxidation (M)",
                             "FragmentMassTolerance": value, "FragmentMassToleranceUnit": unit})
        command = [sys.executable, str(Path(__file__).resolve().parents[1] / "bin/resolve_silac_config_experiment.py"),
                   str(self.openms), str(self.design), "run.raw", str(self.output), str(self.chemistry)]
        if override is not None:
            command.append(f"--comet-fragment-mass-tolerance={override}")
        return subprocess.run(command, capture_output=True, text=True)

    def test_existing_positional_command_and_chemistry_are_preserved(self):
        result = self.run_resolver("0.8", "Da")
        self.assertEqual(result.returncode, 0, result.stderr)
        config = json.loads(self.output.read_text())
        self.assertEqual(config["ffm_labels"], "[][Lys8,Arg10]")
        self.assertEqual(config["binary_modifications"], [0, 1, 1])
        self.assertEqual(config["variable_modifications"], [
            "Oxidation (M)", "Label:13C(6)15N(2) (K)", "Label:13C(6)15N(4) (R)",
        ])
        self.assertEqual(config["comet_fragment_mass_tolerance"]["source"], "metadata")

    def test_optional_named_override_is_recorded(self):
        result = self.run_resolver("20", "ppm", "0.02")
        self.assertEqual(result.returncode, 0, result.stderr)
        config = json.loads(self.output.read_text())
        self.assertEqual(config["openms_design"]["FragmentMassTolerance"], "20")
        self.assertEqual(config["comet_fragment_mass_tolerance"], {
            "value": 0.02, "unit": "Da", "source": "explicit_override",
            "metadata_value": "20", "metadata_unit": "ppm",
        })

    def test_invalid_settings_fail_before_writing_output(self):
        for value, unit, override in (("20", "ppm", None), ("0.8", "Da", "0")):
            with self.subTest(override=override):
                result = self.run_resolver(value, unit, override)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("run.raw:", result.stderr)
                self.assertFalse(self.output.exists())


class ExperimentalDesignLookupTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.design = Path(self.directory.name) / "design.tsv"

    def write_design(self, *filenames):
        rows = ["Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample"]
        for group, filename in enumerate(filenames, start=1):
            for label in (1, 2):
                rows.append(f"{group}\t1\t{filename}\t{label}\t{label}")
        rows.extend(["", "Sample\tMSstats_Condition", "1\tlight", "2\theavy"])
        self.design.write_text("\n".join(rows) + "\n")

    def test_raw_input_matches_converted_design_with_all_channels(self):
        self.write_design("design/run.mzML", "design/other.mzML")
        for suffix in ("raw", "RAW", "RaW"):
            with self.subTest(suffix=suffix):
                rows = read_experimental_design(self.design, f"inputs/run.{suffix}")
                self.assertEqual([row["Label"] for row in rows], ["1", "2"])
                self.assertTrue(all(row["Spectra_Filepath"] == "design/run.mzML" for row in rows))

    def test_original_raw_and_mzml_designs_still_match(self):
        for filename in ("run.raw", "run.RAW", "run.mzML", "run.MZML"):
            with self.subTest(filename=filename):
                self.write_design(filename)
                self.assertEqual(len(read_experimental_design(self.design, filename)), 2)

    def test_exact_raw_match_takes_precedence_over_converted_name(self):
        self.write_design("run.raw", "run.mzML")
        rows = read_experimental_design(self.design, "run.raw")
        self.assertEqual({row["Fraction_Group"] for row in rows}, {"1"})

    def test_missing_run_and_unrelated_extensions_fail(self):
        self.write_design("run.mzML")
        for filename in ("other.raw", "run.mzXML"):
            with self.subTest(filename=filename):
                with self.assertRaisesRegex(ValueError, "No matching row"):
                    read_experimental_design(self.design, filename)


if __name__ == "__main__":
    unittest.main()
