"""Scientific interpretation and graceful degradation for completed-run analysis."""

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import analyze_completed_run as analysis
import test_export_silac_ratios as fixtures


class ScientificAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.RatioExportTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.directory.cleanup)
        self.peptides, self.proteins = self.fixture.export()
        self.groups, _ = analysis.exporter.read_design(self.fixture.design)
        self.labels = ["[]", "[Lys4,Arg6]", "[Lys8,Arg10]"]

    def test_non_light_reference_and_derived_comparisons(self):
        pairs = list(analysis.pair_ratios(self.peptides, self.groups, self.labels, 2))
        raw = [row for row in pairs if row["normalization"] == "raw"]
        self.assertEqual([row["log2_ratio"] for row in raw], [1, 2, 1])
        self.assertEqual([row["estimate_kind"] for row in raw],
                         ["inverse_of_stored", "derived_from_stored_ratios", "stored"])
        self.assertEqual(len(pairs), 6)

    def test_missing_source_does_not_create_derived_ratio(self):
        peptides = [row | {"ratio_status": "not_reported", "log2_ratio": None}
                    if row["numerator_channel"] == 1 else row for row in self.peptides]
        pairs = list(analysis.pair_ratios(peptides, self.groups, self.labels, 2))
        derived = [row for row in pairs if row["estimate_kind"] == "derived_from_stored_ratios"]
        self.assertTrue(all(row["log2_ratio"] is None for row in derived))

    def test_published_pxd039580_mixtures_and_protein_expectations(self):
        profiles = analysis.load_profiles(ROOT / "conf/analysis_profiles.json")
        profile = profiles["PXD039580"]
        groups = {1: {"runs": {(1, "20220603_QE1_TB_A3R1_DDA1")}}}
        pairs = [{"entity": "P1", "sequence": "AKR", "fraction_group": 1, "runs": "mix1",
                  "normalization": "raw", "numerator_label": self.labels[2],
                  "denominator_label": self.labels[0], "estimate_kind": "stored",
                  "log2_ratio": math.log2(15 / 70)}]
        result = analysis.benchmark(pairs, groups, profile, "protein")[0]
        self.assertAlmostEqual(result["expected_ratio"], 15 / 70)
        self.assertAlmostEqual(result["median_log2_error"], 0)
        self.assertEqual(analysis.benchmark(pairs, groups, profiles["PXD003327"], "protein"), [])

    def test_missing_planned_replicate_reduces_completeness_not_n(self):
        groups = {i: {"runs": {(1, f"run{i}")}} for i in (1, 2)}
        profile = {"runs": {f"run{i}": {"replicate_set": "mix", "replicate": str(i)} for i in (1, 2, 3)}}
        pairs = [{"fraction_group": i, "entity": "P1", "normalization": "raw",
                  "numerator_label": "H", "denominator_label": "L", "log2_ratio": value}
                 for i, value in ((1, 0), (2, 1))]
        precision, correlations = analysis.reproducibility(pairs, groups, profile, "protein")
        self.assertEqual(precision[0]["replicates_expected"], 3)
        self.assertEqual(precision[0]["replicates_reported"], 2)
        self.assertFalse(precision[0]["complete"])
        self.assertAlmostEqual(precision[0]["cv_percent"], math.sqrt(.5) / 1.5 * 100)
        self.assertEqual(len(correlations), 3)
        self.assertIsNone(correlations[0]["pearson_log2"])

    def test_label_swaps_and_technical_repeats_collapse_before_biological_test(self):
        groups, profile, pairs = {}, {"runs": {}, "contrasts": [
            {"name": "drug/control", "numerator": "drug", "denominator": "control"}]}, []
        for index, (bio, effect, swapped) in enumerate((("a", .9, False), ("a", 1.1, False),
                                                       ("b", 1.2, True), ("c", .8, False)), 1):
            groups[index] = {"runs": {(1, f"run{index}")}}
            profile["runs"][f"run{index}"] = {
                "biological_replicate": bio,
                "conditions": {"H": "control" if swapped else "drug", "L": "drug" if swapped else "control"}}
            pairs.append({"fraction_group": index, "normalization": "raw", "entity": "P1",
                          "numerator_label": "H", "denominator_label": "L", "log2_ratio": -effect if swapped else effect})
        result = analysis.biological_effects(pairs, groups, profile)[0]
        self.assertEqual(result["biological_replicates_reported"], 3)
        self.assertAlmostEqual(result["mean_log2_fold_change"], 1)
        self.assertAlmostEqual(result["sd_log2"], .2)
        self.assertIn(result["test_status"], ("tested", "scipy_unavailable"))
        if result["p_value"] is not None:
            self.assertEqual(result["q_value"], result["p_value"])

    def test_bh_adjustment_and_undefined_correlations(self):
        rows = [{"p_value": p, "q_value": None} for p in (.01, .04, .03, None)]
        analysis.adjust_bh(rows)
        self.assertEqual([row["q_value"] for row in rows], [.03, .04, .04, None])
        self.assertIsNone(analysis.pearson([1, 1, 1], [1, 2, 3]))
        self.assertIsNone(analysis.pearson([1, 2], [1, 2]))
        self.assertAlmostEqual(analysis.pearson([1, 2, 3], [2, 4, 6]), 1)

    def test_linear_median_is_not_exponentiated_log_median_for_even_n(self):
        self.assertEqual(analysis.median_linear_from_logs([0, 2]), 2.5)
        self.assertIsNone(analysis.median_linear_from_logs([]))
        self.assertIsNone(analysis.median_linear_from_logs([2000]))

    def test_ambiguous_design_is_rejected(self):
        other = self.fixture.folder / "other_design.tsv"
        other.write_text(self.fixture.design.read_text().replace("treatment", "unrelated"))
        with self.assertRaisesRegex(ValueError, "Multiple different designs"):
            analysis.find_export(self.fixture.mztab, self.fixture.consensus, [[self.fixture.design, other]])

    def test_empty_and_partial_result_directories_still_write_reports(self):
        outdir = self.fixture.folder / "analysis"
        self.fixture.consensus.unlink()
        with patch("builtins.print"):
            status = analysis.main([str(self.fixture.folder), "--outdir", str(outdir), "--no-plots", "--strict"])
        self.assertEqual(status, 2)
        metrics = json.loads((outdir / "metrics.json").read_text())
        self.assertEqual(metrics["bundles"][0]["status"], "partial")
        self.assertIsNotNone(metrics["bundles"][0]["identification"])
        self.assertIn("consensusXML absent", metrics["bundles"][0]["diagnostics"][0])
        self.assertTrue((outdir / "report.html").exists())
        with tempfile.TemporaryDirectory() as empty:
            with patch("builtins.print"):
                self.assertEqual(analysis.main([empty, "--no-plots"]), 0)
            self.assertTrue(json.loads((Path(empty) / "analysis/metrics.json").read_text())["diagnostics"])

    def test_recursive_cohorts_choose_local_design_and_escape_html(self):
        outdir = self.fixture.folder / "analysis"
        local = self.fixture.folder / "cohort"
        local.mkdir()
        for path in (self.fixture.mztab, self.fixture.consensus, self.fixture.design):
            (local / ("experimental_design.tsv" if path == self.fixture.design else path.name)).write_bytes(path.read_bytes())
        with patch("builtins.print"):
            self.assertEqual(analysis.main([str(local), "--outdir", str(outdir), "--no-plots"]), 0)
        metrics = json.loads((outdir / "metrics.json").read_text())
        self.assertEqual(metrics["bundles"][0]["status"], "analyzed")
        metrics["bundles"][0]["answers"] = ["<script>alert(1)</script>"]
        self.assertNotIn("<script>", analysis.render_report(metrics, []))

    @unittest.skipUnless((ROOT / "results/PXD003327_duplex/experiment.consensusXML").exists(),
                         "Local duplex result is not available")
    def test_new_analysis_reproduces_saved_benchmark_medians(self):
        baseline = json.loads((ROOT / "tests/fixtures/pxd003327_duplex_baseline.json").read_text())
        with tempfile.TemporaryDirectory() as output:
            with patch("builtins.print"):
                status = analysis.main([str(ROOT / "results/PXD003327_duplex"), "--outdir", output, "--no-plots"])
            self.assertEqual(status, 0)
            bundle = json.loads((Path(output) / "metrics.json").read_text())["bundles"][0]
        self.assertEqual(bundle["status"], "analyzed")
        self.assertEqual(bundle["diagnostics"], [])
        for expected in baseline["benchmark"]:
            observed = next(row for row in bundle["benchmark"] if row["runs"] == expected["run"]
                            and row["peptide_class"] == expected["class"] + "_only")
            self.assertEqual(observed["reported"], expected["n"])
            self.assertAlmostEqual(observed["median_ratio"], expected["median"], places=12)


if __name__ == "__main__":
    unittest.main()
