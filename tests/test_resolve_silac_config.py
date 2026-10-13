"""Design lookup before RAW inputs have been converted to mzML."""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from resolve_silac_config_experiment import read_experimental_design


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
