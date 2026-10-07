"""Small format fixtures plus an optional regression against the local duplex result."""

import csv
from collections import Counter
import importlib.util
import json
import math
import statistics
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("exporter", ROOT / "bin/export_silac_ratios.py")
EXPORTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORTER)
PREFIX = EXPORTER.PEPTIDE_PREFIX


class RatioExportTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.folder = Path(self.directory.name)
        self.mztab = self.folder / "experiment.mzTab"
        self.consensus = self.folder / "experiment.consensusXML"
        self.design = self.folder / "design.tsv"
        lines = ["Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample"]
        for fraction in (2, 1):
            for channel in (1, 2, 3):
                lines.append(f"7\t{fraction}\trun_F{fraction}.mzML\t{channel}\t{channel}")
        lines += ["", "Sample\tMSstats_Condition\tMSstats_BioReplicate",
                  "1\tcontrol\t1", "2\tstandard\t1", "3\ttreatment\t1"]
        self.design.write_text("\n".join(lines) + "\n")
        self.root = ET.Element("consensusXML")
        columns = ET.SubElement(self.root, "mapList")
        for fraction in (1, 2):
            for channel in (1, 2, 3):
                column = ET.SubElement(columns, "map", name=f"run_F{fraction}.mzML")
                self.params(column, {"fraction_group": "7", "fraction": str(fraction),
                                     "channel_id": str(channel - 1)})
        processing = ET.SubElement(self.root, "dataProcessing")
        ET.SubElement(processing, "software", name="MS1LabeledWorkflow")
        self.params(processing, {
            "parameter: labels": "[][Lys4,Arg6][Lys8,Arg10]",
            "parameter: ratios:reference_channel": "2",
            "parameter: ratios:normalize": "true",
        })
        protein = ET.SubElement(ET.SubElement(self.root, "IdentificationRun"), "ProteinIdentification")
        ET.SubElement(protein, "ProteinHit", id="PH_A", accession="P1")
        ET.SubElement(protein, "ProteinHit", id="PH_B", accession="P2")
        self.params(protein, {
            "indistinguishable_proteins_0": "0.001,PH_B,PH_A",
            "indistinguishable_proteins_0_fraction_group_level_ratio": "[0.5,1,2]",
            "indistinguishable_proteins_0_fraction_group_level_ratio_normalized": "[0.25,1,1.5]",
            "indistinguishable_proteins_0_fraction_group_level_ratio_fraction_group": "[7,7,7]",
            "indistinguishable_proteins_0_fraction_group_level_ratio_label": "[1,2,3]",
            "indistinguishable_proteins_0_fraction_group_level_ratio_count": "[3,3,3]",
            # Same leader, different full membership; no stored ratios.
            "indistinguishable_proteins_1": "0.001,PH_A",
        })
        self.peptide = {
            "sequence": "AKR", "modifications": "null", "charge": "2",
            "opt_global_cv_MS:1000889_peptidoform_sequence": "AKR",
            "opt_global_cv_MS:1002217_decoy_peptide": "0",
            PREFIX: "[0.5,1,2]", PREFIX + "_normalized": "[0.25,1,1.5]",
            PREFIX + "_fraction_group": "[7,7,7]", PREFIX + "_channel": "[1,2,3]",
            PREFIX + "_count": "[4,4,4]",
        }

    @staticmethod
    def params(node, values):
        for name, value in values.items():
            ET.SubElement(node, "UserParam", name=name, value=value)

    def export(self, rows=None):
        ET.ElementTree(self.root).write(self.consensus)
        with self.mztab.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["PEH", *self.peptide])
            for row in rows if rows is not None else [self.peptide]:
                writer.writerow(["PEP", *(row[field] for field in self.peptide)])
        return EXPORTER.export_ratios(self.mztab, self.consensus, self.design)

    def test_stored_ratios_full_membership_fractions_and_reference(self):
        peptides, proteins = self.export()
        raw = [row for row in peptides if row["normalization"] == "raw"]
        self.assertEqual([row["ratio"] for row in raw], [0.5, 2.0])
        self.assertEqual([row["log2_ratio"] for row in raw], [-1.0, 1.0])
        self.assertTrue(all(row["runs"] == "run_F1;run_F2" for row in peptides))
        self.assertTrue(all(row["reference_channel"] == 2 for row in peptides))
        self.assertTrue(all(row["reference_label_definition"] == "[Lys4,Arg6]" for row in peptides))
        self.assertEqual(raw[0]["numerator_condition"], "control")
        self.assertEqual({row["protein_group"] for row in proteins}, {"P1", "P1;P2"})
        self.assertEqual([row["ratio_count"] for row in proteins if row["ratio_status"] == "reported"], [3]*4)
        normalized = [row for row in peptides if row["normalization"] == "normalized"]
        self.assertEqual([row["ratio"] for row in normalized], [0.25, 1.5])

    def test_charge_duplicates_deduplicated_and_decoys_excluded(self):
        peptides, _ = self.export([self.peptide, self.peptide | {"charge": "3"},
                                   self.peptide | {"sequence": "DECOYK", "opt_global_cv_MS:1002217_decoy_peptide": "1"}])
        self.assertEqual(len(peptides), 4)
        self.assertEqual({row["sequence"] for row in peptides}, {"AKR"})

    def test_modifications_remain_distinct(self):
        peptides, _ = self.export([self.peptide, self.peptide | {
            "modifications": "1-UNIMOD:35", "opt_global_cv_MS:1000889_peptidoform_sequence": "A(Oxidation)KR"}])
        self.assertEqual(len(peptides), 8)

    def test_missing_quantities_preserved(self):
        empty = self.peptide | {PREFIX: "null", PREFIX + "_normalized": "null"}
        peptides, proteins = self.export([empty])
        self.assertTrue(all(row["ratio"] is None and row["ratio_status"] == "not_reported" for row in peptides))
        self.assertTrue(all(row["ratio_status"] == "not_reported" for row in proteins if row["protein_group"] == "P1"))

    def test_invalid_values_have_no_log2(self):
        peptides, _ = self.export([self.peptide | {PREFIX: "[0,1,nan]"}])
        raw = [row for row in peptides if row["normalization"] == "raw"]
        self.assertEqual([row["ratio_status"] for row in raw], ["non_positive", "non_finite"])
        self.assertTrue(all(row["log2_ratio"] is None for row in raw))

    def test_conflicting_duplicates_fail(self):
        with self.assertRaisesRegex(ValueError, "Conflicting repeated ratio"):
            self.export([self.peptide, self.peptide | {PREFIX: "[0.75,1,2]"}])

    def test_misaligned_arrays_fail(self):
        with self.assertRaisesRegex(ValueError, "different lengths"):
            self.export([self.peptide | {PREFIX + "_count": "[4,4]"}])

    def test_missing_design_channel_fails(self):
        with self.assertRaisesRegex(ValueError, "absent from the design"):
            self.export([self.peptide | {PREFIX + "_channel": "[1,2,4]"}])

    def test_mismatched_design_fails(self):
        self.design.write_text(self.design.read_text().replace("run_F1", "wrong_run"))
        with self.assertRaisesRegex(ValueError, "does not match consensusXML"):
            self.export()


class DuplexRegressionTests(unittest.TestCase):
    @unittest.skipUnless((ROOT / "results/PXD003327_duplex/experiment.consensusXML").exists(),
                         "Local duplex result is not available")
    def test_saved_duplex_baseline(self):
        baseline = json.loads((ROOT / "tests/fixtures/pxd003327_duplex_baseline.json").read_text())
        folder = ROOT / "results/PXD003327_duplex"
        with (folder / "experiment.mzTab").open() as handle:
            counts = Counter(row[0] for row in csv.reader(handle, delimiter="\t")
                             if row and row[0] in ("PRT", "PEP", "PSM"))
        self.assertEqual(dict(counts), baseline["mztab_rows"])
        peptides, proteins = EXPORTER.export_ratios(
            folder / "experiment.mzTab", folder / "experiment.consensusXML",
            folder / "PXD003327_duplex.sdrf_openms_design.tsv")
        peptides = [row for row in peptides if row["ratio_status"] == "reported"]
        proteins = [row for row in proteins if row["ratio_status"] == "reported"]
        self.assertEqual(len(peptides), baseline["peptide_ratio_entries"])
        self.assertEqual(len(proteins), baseline["protein_ratio_entries"])
        self.assertEqual(len({row["peptidoform"] for row in peptides}), baseline["unique_peptides_with_ratio"])
        self.assertEqual(len({row["protein_group"] for row in proteins}), baseline["protein_groups_with_ratio"])
        for expected in baseline["benchmark"]:
            other = "R" if expected["class"] == "K" else "K"
            selected = [row["ratio"] for row in peptides
                        if row["fraction_group"] == expected["fraction_group"]
                        and row["sequence"].count(expected["class"]) == 1
                        and other not in row["sequence"]]
            self.assertEqual(len(selected), expected["n"])
            self.assertTrue(math.isclose(statistics.median(selected), expected["median"], rel_tol=1e-12))


if __name__ == "__main__":
    unittest.main()
