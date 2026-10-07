"""QC semantics on small format fixtures and the optional local duplex baseline."""

import csv
import json
import math
import statistics
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import silac_qc as qc
import export_silac_ratios as exporter
import test_export_silac_ratios as fixtures


class QCStatisticsTests(unittest.TestCase):
    def test_missing_values_stay_in_denominator_and_quartiles_are_inclusive(self):
        rows = [{"ratio_status": "reported", "ratio": str(value), "ratio_count": str(count)}
                for value, count in ((1, 1), (2, 3), (4, 5))]
        rows.append({"ratio_status": "not_reported", "ratio": "null", "ratio_count": "null"})
        summary = qc.distribution(rows)
        self.assertEqual(summary["total_entities"], 4)
        self.assertEqual(summary["reported"], 3)
        self.assertEqual(summary["coverage"], 0.75)
        self.assertEqual(summary["median_ratio"], 2)
        self.assertEqual(summary["median_log2_ratio"], 1)
        self.assertEqual(summary["iqr_log2"], 1)
        self.assertEqual(summary["ratio_count_median"], 3)

    def test_empty_and_singleton_distributions(self):
        self.assertIsNone(qc.distribution([])["median_ratio"])
        singleton = qc.distribution([{"ratio_status": "reported", "ratio": "4", "ratio_count": "null"}])
        self.assertEqual(singleton["median_log2_ratio"], 2)
        self.assertEqual(singleton["iqr_log2"], 0)
        self.assertIsNone(singleton["ratio_count_median"])

    def test_nonpositive_and_nonfinite_excluded_from_distribution(self):
        rows = [{"ratio_status": status, "ratio": "null", "ratio_count": "null"}
                for status in ("non_positive", "non_finite", "not_reported")]
        summary = qc.distribution(rows)
        self.assertEqual(summary["coverage"], 0)
        self.assertEqual(summary["non_positive"], 1)
        self.assertIsNone(summary["median_log2_ratio"])


class QCInputTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.RatioExportTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.directory.cleanup)
        peptides, proteins = self.fixture.export()
        self.peptides = self.fixture.folder / "peptide_ratios.tsv"
        self.proteins = self.fixture.folder / "protein_ratios.tsv"
        exporter.write_table(self.peptides, exporter.PEPTIDE_FIELDS, peptides)
        exporter.write_table(self.proteins, exporter.PROTEIN_FIELDS, proteins)

    def analyze(self, config=None):
        return qc.analyze(self.peptides, self.proteins, self.fixture.mztab,
                          self.fixture.consensus, self.fixture.design, config)

    def test_triplex_reference_and_normalization_are_separate(self):
        metrics, _, _ = self.analyze()
        self.assertEqual(len(metrics["comparisons"]), 8)
        selected = {row["numerator_channel"]: row for row in metrics["comparisons"]
                    if row["level"] == "peptide" and row["normalization"] == "raw"}
        self.assertEqual(selected["1"]["median_ratio"], 0.5)
        self.assertEqual(selected["3"]["median_ratio"], 2)
        self.assertEqual(selected["3"]["reference_channel"], "2")
        self.assertEqual(selected["3"]["runs"], "run_F1;run_F2")
        self.assertEqual(metrics["entities"]["protein_groups_in_export"], 2)
        self.assertEqual(metrics["entities"]["protein_groups_with_raw_ratio"], 1)
        self.assertFalse(metrics["benchmark"])
        json.dumps(metrics, allow_nan=False)

    def test_benchmark_supports_fraction_groups_and_non_light_reference(self):
        config = {"name": "fixture", "peptide_classes": {"one_K_one_R": {"K": 1, "R": 1}},
                  "comparisons": [{"runs": ["run_F2", "run_F1"], "numerator_label_definition": "[]",
                                   "reference_label_definition": "[Lys4,Arg6]", "normalization": "raw",
                                   "expected_ratios": {"one_K_one_R": 0.5}}]}
        metrics, _, _ = self.analyze(config)
        result = metrics["benchmark"][0]
        self.assertEqual(len(metrics["benchmark"]), 1)
        self.assertEqual(result["reported"], 1)
        self.assertEqual(result["median_bias_percent"], 0)
        self.assertEqual(result["median_log2_error"], 0)

    def test_stale_tables_are_rejected(self):
        text = self.peptides.read_text().replace("\traw\t0.5\t", "\traw\t0.75\t")
        self.peptides.write_text(text)
        with self.assertRaisesRegex(ValueError, "differs from its scientific sources"):
            self.analyze()

    def test_report_escapes_metadata(self):
        metrics, _, _ = self.analyze()
        metrics["benchmark_name"] = '<script>alert("unsafe")</script>'
        report = qc.render_report(metrics, [])
        self.assertNotIn('<script>alert("unsafe")</script>', report)
        self.assertIn("&lt;script&gt;", report)
        self.assertIn("all exported experiment-wide entities", report)

    def test_psms_deduplicated_by_spectrum_and_peptide_and_pep_is_not_qvalue(self):
        self.fixture.mztab.write_text(
            "MTD\tms_run[1]-location\tfile://run.mzML\n"
            "MTD\tpsm_search_engine_score[1]\t[MS, MS:1001493, posterior error probability, ]\n"
            "PSH\tspectra_ref\tsequence\tmodifications\tcharge\taccession\tsearch_engine_score[1]\topt_global_cv_MS:1002217_decoy_peptide\n"
            "PSM\tms_run[1]:scan=1\tAKR\tnull\t2\tP1\t0.002\t0\n"
            "PSM\tms_run[1]:scan=1\tAKR\tnull\t2\tP2\t0.002\tnull\n")
        result, _ = qc.read_identifications(self.fixture.mztab)
        self.assertEqual(result["mztab_rows"]["PSM"], 2)
        self.assertEqual(result["unique_resolved_psms"], 1)
        self.assertEqual(result["psm_decoy_annotations"]["unknown_rows"], 1)
        self.assertFalse(result["psm_score"]["is_qvalue"])
        self.assertIsNone(result["psm_score"]["max_reported_qvalue"])


class QCDuplexRegressionTests(unittest.TestCase):
    @unittest.skipUnless((ROOT / "results/PXD003327_duplex/tables/peptide_ratios.tsv").exists(),
                         "Local duplex ratio tables are not available")
    def test_local_duplex_matches_preserved_baseline(self):
        folder = ROOT / "results/PXD003327_duplex"
        baseline = json.loads((ROOT / "tests/fixtures/pxd003327_duplex_baseline.json").read_text())
        metrics, _, _ = qc.analyze(folder / "tables/peptide_ratios.tsv", folder / "tables/protein_ratios.tsv",
                                   folder / "experiment.mzTab", folder / "experiment.consensusXML",
                                   folder / "PXD003327_duplex.sdrf_openms_design.tsv",
                                   qc.read_benchmark(ROOT / "assets/benchmarks/pxd003327.json"))
        self.assertEqual(metrics["identification"]["mztab_rows"], baseline["mztab_rows"])
        self.assertEqual(metrics["entities"]["peptides_with_raw_ratio"], 701)
        self.assertEqual(metrics["entities"]["protein_groups_with_raw_ratio"], 117)
        self.assertEqual(sum(row["reported"] for row in metrics["comparisons"] if row["level"] == "peptide"), 954)
        self.assertEqual(sum(row["reported"] for row in metrics["comparisons"] if row["level"] == "protein"), 182)
        for expected in baseline["benchmark"]:
            observed = next(row for row in metrics["benchmark"]
                            if row["runs"] == expected["run"] and row["peptide_class"] == expected["class"] + "_only")
            self.assertEqual(observed["reported"], expected["n"])
            self.assertTrue(math.isclose(observed["median_ratio"], expected["median"], rel_tol=1e-12))


if __name__ == "__main__":
    unittest.main()
