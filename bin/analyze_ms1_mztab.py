#!/usr/bin/env python3
"""Summarize MS1LabeledWorkflow peptide ratios from mzTab.

Usage:
  python3 bin/analyze_ms1_mztab.py experiment.mzTab experimental_design.tsv
  python3 bin/analyze_ms1_mztab.py experiment.mzTab experimental_design.tsv --benchmark pxd003327

Ratios come from OpenMS peptide_ratio fields (channel / light), not from
the assay abundance columns. Repeated mzTab peptide rows are counted once.
"""

import argparse
import ast
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED = {
    "Chris_Ecoli_4-1": {"K": {"H/L": 0.25}, "R": {"H/L": 4.0}},
    "Chris_Ecoli_1-1": {"K": {"H/L": 1.0}, "R": {"H/L": 1.0}},
    "Chris_Ecoli_1-2-4": {
        "K": {"M/L": 2.0, "H/L": 4.0},
        "R": {"M/L": 0.5, "H/L": 0.25},
    },
}


def read_design(path):
    with path.open(newline="") as handle:
        tables = ([], [])
        section = 0
        for line in handle:
            if not line.strip():
                section = 1
                continue
            tables[section].append(line)

    runs = defaultdict(set)
    channels = {}
    samples = {
        row["Sample"]: row["MSstats_Condition"].lower()
        for row in csv.DictReader(tables[1], delimiter="\t")
    }
    for row in csv.DictReader(tables[0], delimiter="\t"):
        group = row["Fraction_Group"]
        channel = int(row["Label"])
        runs[group].add(Path(row["Spectra_Filepath"]).stem)
        condition = samples[row["Sample"]]
        name = next((name for name in ("light", "medium", "heavy") if name in condition), None)
        if name is None:
            raise ValueError(f"Unknown SILAC condition for sample {row['Sample']}: {condition}")
        channels[(group, channel)] = name
    if not runs:
        raise ValueError(f"No runs found in design: {path}")
    return {group: sorted(names) for group, names in runs.items()}, channels


def peptide_class(sequence):
    if sequence.count("K") == 1 and "R" not in sequence:
        return "K"
    if sequence.count("R") == 1 and "K" not in sequence:
        return "R"
    return None


def read_ratios(path):
    counts = Counter()
    ratios = defaultdict(set)
    header = None
    fields = (
        "opt_global_MS1Label:peptide_ratio",
        "opt_global_MS1Label:peptide_ratio_channel",
        "opt_global_MS1Label:peptide_ratio_fraction_group",
    )

    with path.open(newline="") as handle:
        for record in csv.reader(handle, delimiter="\t"):
            if not record:
                continue
            if record[0] == "PEH":
                header = record
                missing = set(fields) - set(header)
                if missing:
                    raise ValueError(f"Missing mzTab peptide fields: {sorted(missing)}")
            elif record[0] in ("PRT", "PEP", "PSM"):
                counts[record[0]] += 1
                if record[0] != "PEP":
                    continue
                if header is None:
                    raise ValueError("PEP row found before PEH header")
                row = dict(zip(header, record))
                if row.get("opt_global_cv_MS:1002217_decoy_peptide") == "1":
                    continue
                kind = peptide_class(row["sequence"])
                if kind is None or row[fields[0]] == "null":
                    continue
                values = [ast.literal_eval(row[field]) for field in fields]
                if len({len(value) for value in values}) != 1:
                    raise ValueError("Peptide ratio arrays have different lengths")
                for ratio, channel, group in zip(*values):
                    ratio = float(ratio)
                    channel = int(channel)
                    if channel > 1 and math.isfinite(ratio) and ratio > 0:
                        key = (row["sequence"], row["modifications"])
                        ratios[(str(group), channel, kind)].add((key, ratio))

    if header is None:
        raise ValueError(f"No peptide section found in mzTab: {path}")
    if not ratios:
        raise ValueError(f"No positive SILAC peptide ratios found in mzTab: {path}")
    return counts, ratios


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mztab", type=Path, help="MS1LabeledWorkflow experiment.mzTab")
    parser.add_argument("design", type=Path, help="matching OpenMS experimental design TSV")
    parser.add_argument("--benchmark", choices=["pxd003327"], help="show known PXD003327 ratios")
    args = parser.parse_args()

    runs, channels = read_design(args.design)
    counts, ratios = read_ratios(args.mztab)
    print("Records: " + ", ".join(f"{tag}={counts[tag]}" for tag in ("PRT", "PEP", "PSM")))
    print("Group\tRun\tClass\tRatio\tN\tMedian\tExpected")
    for (group, channel, kind), observations in sorted(ratios.items()):
        names = runs.get(group)
        if not names:
            raise ValueError(f"Fraction group {group} is absent from the design")
        run = ",".join(names)
        reference = channels.get((group, 1))
        numerator = channels.get((group, channel))
        if reference != "light" or numerator not in ("medium", "heavy"):
            raise ValueError(f"Unexpected channel mapping for group {group}: {reference}, {numerator}")
        label = f"{numerator[0].upper()}/L"
        values = [ratio for _, ratio in observations]
        median = statistics.median(values)
        expected = ""
        if args.benchmark:
            if len(names) != 1:
                raise ValueError(f"Benchmark needs one run in fraction group {group}")
            target = EXPECTED.get(run, {}).get(kind, {}).get(label)
            expected = f"{target:.3f}" if target is not None else ""
        print(f"{group}\t{run}\t{kind}\t{label}\t{len(values)}\t{median:.3f}\t{expected}")


if __name__ == "__main__":
    main()
