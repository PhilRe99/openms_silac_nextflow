#!/usr/bin/env python3
"""Analyze published OpenMS SILAC runs, including incomplete result directories.

Python >=3.10; standard-library analysis, optional Matplotlib plots and SciPy
paired biological tests. See docs/completed_run_analysis.md for definitions.
"""

import argparse
import base64
import csv
import hashlib
import html
import itertools
import importlib.metadata
import json
import math
import re
import statistics
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import export_silac_ratios as exporter
import silac_qc as qc


ROOT = Path(__file__).resolve().parents[1]
ERRORS = (ValueError, KeyError, TypeError, IndexError, OSError, ET.ParseError)
BENCHMARK_FIELDS = ["level", "fraction_group", "runs", "peptide_class", "normalization",
                    "numerator_label", "denominator_label", "estimate_kind", "expected_ratio",
                    "total_entities", "reported", "coverage", "median_ratio", "median_log2_error",
                    "median_bias_percent", "median_absolute_log2_error", "iqr_log2"]
PRECISION_FIELDS = ["level", "replicate_set", "replicate_type", "normalization", "numerator_label",
                    "denominator_label", "entity", "replicates_expected", "replicates_reported",
                    "complete", "cv_percent", "sd_log2"]
CORRELATION_FIELDS = ["level", "replicate_set", "normalization", "numerator_label",
                      "denominator_label", "replicate_a", "replicate_b", "shared_entities",
                      "pearson_log2", "median_absolute_log2_difference"]
EFFECT_FIELDS = ["contrast", "protein_group", "normalization", "biological_replicates_expected",
                 "biological_replicates_reported", "mean_log2_fold_change", "sd_log2",
                 "p_value", "q_value", "test_status"]


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path.resolve()), "sha256": digest.hexdigest()}


def table(path, rows, fields):
    qc.write_summary(path, rows, fields)


def text_rows(rows):
    return [{key: "null" if value is None else str(value) for key, value in row.items()} for row in rows]


def validate_profile(profile):
    if not isinstance(profile.get("runs", {}), dict):
        raise ValueError("Profile runs must be a mapping keyed by exact run stems")
    for name, counts in profile.get("peptide_classes", {}).items():
        if any(not re.fullmatch("[A-Z]", residue) or type(count) is not int or count < 0
               for residue, count in counts.items()):
            raise ValueError(f"Invalid peptide class: {name}")
    for run, metadata in profile.get("runs", {}).items():
        if not isinstance(metadata, dict):
            raise ValueError(f"Run metadata must be an object: {run}")
        if metadata.get("replicate_type", "technical") not in ("technical", "biological"):
            raise ValueError(f"Invalid replicate_type: {run}")
        if metadata.get("replicate_set") and not metadata.get("replicate"):
            raise ValueError(f"replicate_set requires replicate identity: {run}")
        for kind, amounts in metadata.get("amounts", {}).items():
            if kind not in profile.get("peptide_classes", {}):
                raise ValueError(f"Amounts reference undefined peptide class: {kind}")
            if not amounts or any(not math.isfinite(float(value)) or float(value) <= 0
                                  for value in amounts.values()):
                raise ValueError(f"Mixture amounts must be finite and positive: {run}/{kind}")
    names = []
    for contrast in profile.get("contrasts", []):
        names.append(contrast["name"])
        if not contrast["numerator"] or not contrast["denominator"] or contrast["numerator"] == contrast["denominator"]:
            raise ValueError("Contrasts require distinct numerator and denominator conditions")
    if len(names) != len(set(names)):
        raise ValueError("Contrast names must be unique")
    return profile


def load_profiles(path):
    data = json.loads(path.read_text())
    if data.get("schema_version") != 1:
        raise ValueError("Unsupported analysis profile schema_version")
    return {name: validate_profile(profile) for name, profile in data["datasets"].items()}


def choose_profile(profiles, dataset, mztab, runs):
    if dataset == "none":
        return "unknown", {}
    if dataset != "auto":
        if dataset not in profiles:
            raise ValueError(f"Dataset {dataset!r} is absent from the profiles")
        return dataset, profiles[dataset]
    # Exact known run stems provide a fallback when result folders are renamed.
    named = [name for name in profiles if name.lower() in str(mztab).lower()]
    matched = [name for name, profile in profiles.items()
               if runs and runs.issubset(profile.get("runs", {}))
               and profile.get("detect_by_runs", False)]
    candidates = set(named or matched)
    if len(candidates) == 1:
        name = candidates.pop()
        return name, profiles[name]
    return "unknown", {}


def group_metadata(group, profile):
    metadata = [profile.get("runs", {}).get(run) for _, run in sorted(group["runs"])]
    if not metadata or any(item is None for item in metadata):
        return None
    if any(item != metadata[0] for item in metadata[1:]):
        raise ValueError("Fractions in one fraction group have inconsistent scientific metadata")
    return metadata[0]


def pair_ratios(rows, groups, labels, reference):
    """All channel pairs in log space; keep derived ratios explicitly marked.

    Dividing two aggregated protein ratios is a derived estimate, not OpenMS
    quantification with a different reference. Reference/reference is an identity.
    """
    entities = defaultdict(dict)
    for row in rows:
        entity = row.get("protein_group") or (row["sequence"], row["modifications"])
        key = (int(row["fraction_group"]), row["normalization"], entity)
        entities[key][int(row["numerator_channel"])] = row
    for (group_id, normalization, entity), channels in sorted(entities.items()):
        group = groups[group_id]
        for denominator, numerator in itertools.combinations(sorted(group["samples"]), 2):
            logs = {reference: 0.0}
            logs.update({channel: float(row["log2_ratio"]) for channel, row in channels.items()
                         if row["ratio_status"] == "reported"})
            value = logs[numerator] - logs[denominator] if numerator in logs and denominator in logs else None
            kind = ("stored" if denominator == reference else
                    "inverse_of_stored" if numerator == reference else "derived_from_stored_ratios")
            sample = next(iter(channels.values()))
            yield {"fraction_group": group_id, "normalization": normalization,
                   "entity": entity, "sequence": sample.get("sequence", ""),
                   "runs": sample["runs"], "numerator_label": labels[numerator - 1],
                   "denominator_label": labels[denominator - 1], "estimate_kind": kind,
                   "log2_ratio": value}


def benchmark(pairs, groups, profile, level):
    selected = defaultdict(list)
    metadata_by_group = {gid: group_metadata(group, profile) for gid, group in groups.items()}
    for row in pairs:
        metadata = metadata_by_group[row["fraction_group"]]
        if metadata is None:
            continue
        for kind, counts in profile.get("peptide_classes", {}).items():
            # Residue-specific mixtures do not have a universal protein expectation.
            if level == "protein" and counts:
                continue
            if level == "peptide" and any(row["sequence"].count(residue) != count for residue, count in counts.items()):
                continue
            amounts = metadata.get("amounts", {}).get(kind, {})
            num, den = row["numerator_label"], row["denominator_label"]
            if num not in amounts or den not in amounts:
                continue
            expected_log = math.log2(float(amounts[num])) - math.log2(float(amounts[den]))
            key = (row["fraction_group"], row["runs"], kind, row["normalization"], num, den,
                   row["estimate_kind"], expected_log)
            selected[key].append(row["log2_ratio"])
    summaries = []
    for key, values in sorted(selected.items()):
        gid, runs, kind, normalization, num, den, estimate, expected_log = key
        logs = sorted(value for value in values if value is not None)
        errors = [value - expected_log for value in logs]
        median_ratio = median_linear_from_logs(logs)
        relative_median = median_linear_from_logs(errors)
        summaries.append(dict(zip(BENCHMARK_FIELDS[:9],
                                  [level, gid, runs, kind, normalization, num, den, estimate, safe_exp2(expected_log)])) |
                         {"total_entities": len(values), "reported": len(logs), "coverage": len(logs) / len(values),
                          "median_ratio": median_ratio,
                          "median_log2_error": statistics.median(errors) if errors else None,
                          "median_bias_percent": (relative_median - 1) * 100 if relative_median is not None else None,
                          "median_absolute_log2_error": statistics.median(map(abs, errors)) if errors else None,
                          "iqr_log2": qc.quantile(logs, .75) - qc.quantile(logs, .25) if logs else None})
    return summaries


def safe_exp2(value):
    if value is None:
        return None
    try:
        result = 2.0 ** value
    except OverflowError:
        return None
    return result if math.isfinite(result) and result > 0 else None


def median_linear_from_logs(sorted_logs):
    if not sorted_logs:
        return None
    size = len(sorted_logs)
    low = safe_exp2(sorted_logs[(size - 1) // 2])
    high = safe_exp2(sorted_logs[size // 2])
    return low + (high - low) / 2 if low is not None and high is not None else None


def pearson(x, y):
    if len(x) < 3:
        return None
    mx, my = statistics.mean(x), statistics.mean(y)
    dx, dy = [value - mx for value in x], [value - my for value in y]
    denominator = math.sqrt(sum(value * value for value in dx) * sum(value * value for value in dy))
    return max(-1, min(1, sum(a * b for a, b in zip(dx, dy)) / denominator)) if denominator else None


def reproducibility(pairs, groups, profile, level):
    sets = defaultdict(lambda: defaultdict(dict))
    for row in pairs:
        metadata = group_metadata(groups[row["fraction_group"]], profile)
        if not metadata or not metadata.get("replicate_set"):
            continue
        key = (metadata["replicate_set"], metadata.get("replicate_type", "technical"),
               row["normalization"], row["numerator_label"], row["denominator_label"])
        replicate = str(metadata["replicate"])
        entity = row["entity"]
        if entity in sets[key][replicate]:
            raise ValueError("Multiple fraction groups map to one replicate; combine fractions in the design or use distinct replicate IDs")
        sets[key][replicate][entity] = row["log2_ratio"]
    precision, correlations = [], []
    for (name, kind, normalization, num, den), replicates in sorted(sets.items()):
        # Keep planned but absent injections in completeness denominators.
        for metadata in profile.get("runs", {}).values():
            if metadata.get("replicate_set") == name and metadata.get("replicate_type", "technical") == kind:
                replicates.setdefault(str(metadata["replicate"]), {})
        common = {"level": level, "replicate_set": name, "replicate_type": kind,
                  "normalization": normalization, "numerator_label": num, "denominator_label": den}
        for entity in sorted(set().union(*(set(values) for values in replicates.values()))):
            logs = [values[entity] for values in replicates.values() if values.get(entity) is not None]
            # Scaling all ratios by a constant leaves CV unchanged and avoids overflow.
            scaled = [2 ** (value - max(logs)) for value in logs] if logs else []
            precision.append(common | {"entity": ";".join(entity) if isinstance(entity, tuple) else entity,
                                       "replicates_expected": len(replicates), "replicates_reported": len(logs),
                                       "complete": len(logs) == len(replicates),
                                       "cv_percent": statistics.stdev(scaled) / statistics.mean(scaled) * 100 if len(logs) >= 2 else None,
                                       "sd_log2": statistics.stdev(logs) if len(logs) >= 2 else None})
        for a, b in itertools.combinations(sorted(replicates), 2):
            shared = sorted(entity for entity in replicates[a].keys() & replicates[b].keys()
                            if replicates[a][entity] is not None and replicates[b][entity] is not None)
            x, y = [replicates[a][entity] for entity in shared], [replicates[b][entity] for entity in shared]
            correlations.append(common | {"replicate_a": a, "replicate_b": b, "shared_entities": len(shared),
                                          "pearson_log2": pearson(x, y),
                                          "median_absolute_log2_difference": statistics.median(abs(u - v) for u, v in zip(x, y)) if shared else None})
    return precision, correlations


def adjust_bh(rows):
    tested = sorted((row["p_value"], index) for index, row in enumerate(rows) if row["p_value"] is not None)
    previous = 1.0
    for rank in range(len(tested), 0, -1):
        pvalue, index = tested[rank - 1]
        previous = min(previous, pvalue * len(tested) / rank)
        rows[index]["q_value"] = previous


def biological_effects(pairs, groups, profile):
    """Paired within-mixture contrasts only; technical repeats collapse first."""
    try:
        from scipy.stats import ttest_1samp
    except ImportError:
        ttest_1samp = None
    results = []
    for contrast in profile.get("contrasts", []):
        values = defaultdict(lambda: defaultdict(list))
        expected = defaultdict(set)
        for row in pairs:
            metadata = group_metadata(groups[row["fraction_group"]], profile)
            if not metadata or not metadata.get("biological_replicate"):
                continue
            if metadata.get("comparison", "") != contrast.get("comparison", ""):
                continue
            conditions = metadata.get("conditions", {})
            num, den = conditions.get(row["numerator_label"]), conditions.get(row["denominator_label"])
            if (num, den) == (contrast["numerator"], contrast["denominator"]):
                direction = 1
            elif (den, num) == (contrast["numerator"], contrast["denominator"]):
                direction = -1
            else:
                continue
            bio = str(metadata["biological_replicate"])
            expected[row["normalization"]].add(bio)
            key = (row["normalization"], row["entity"])
            values[key]  # Include entities with no estimate.
            if row["log2_ratio"] is not None:
                values[key][bio].append(direction * row["log2_ratio"])
        contrast_rows = []
        planned = {str(metadata["biological_replicate"]) for metadata in profile.get("runs", {}).values()
                   if metadata.get("biological_replicate")
                   and metadata.get("comparison", "") == contrast.get("comparison", "")
                   and {contrast["numerator"], contrast["denominator"]}.issubset(metadata.get("conditions", {}).values())}
        for (normalization, entity), replicates in sorted(values.items()):
            logs = [statistics.median(repeated) for repeated in replicates.values()]
            status, pvalue = "insufficient_biological_replicates", None
            if len(logs) >= 3:
                if statistics.stdev(logs) == 0:
                    status = "zero_variance_not_tested"
                elif ttest_1samp is None:
                    status = "scipy_unavailable"
                else:
                    pvalue = float(ttest_1samp(logs, popmean=0).pvalue)
                    status = "tested" if math.isfinite(pvalue) else "nonfinite_test"
                    if not math.isfinite(pvalue):
                        pvalue = None
            contrast_rows.append({"contrast": contrast["name"], "protein_group": entity,
                                  "normalization": normalization,
                                  "biological_replicates_expected": len(planned | expected[normalization]),
                                  "biological_replicates_reported": len(logs),
                                  "mean_log2_fold_change": statistics.mean(logs) if logs else None,
                                  "sd_log2": statistics.stdev(logs) if len(logs) >= 2 else None,
                                  "p_value": pvalue, "q_value": None, "test_status": status})
        for normalization in {row["normalization"] for row in contrast_rows}:
            adjust_bh([row for row in contrast_rows if row["normalization"] == normalization])
        results.extend(contrast_rows)
    return results


def scientific_plots(outdir, result):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figures = []
    for name, rows, field, xlabel in (
        ("benchmark_error", result["benchmark"], "median_log2_error", "Median log2(observed / expected)"),
        ("replicate_cv", result["precision_summary"], "median_cv_percent_complete", "Median complete-case ratio CV (%)"),
    ):
        selected = [row for row in rows if row[field] is not None]
        if not selected:
            continue
        labels = [f"{row['level']} · {row.get('runs', row.get('replicate_set'))} · "
                  f"{row['numerator_label']}/{row['denominator_label']} · {row['normalization']}"
                  + (f" · {row['peptide_class']}" if "peptide_class" in row else "") for row in selected]
        figure, axis = plt.subplots(figsize=(11, max(3, len(selected) * .32)))
        axis.scatter([row[field] for row in selected], range(len(selected)), color="#287e8a")
        axis.axvline(0, color="#b75337", linestyle="--")
        axis.set_yticks(range(len(selected)), labels, fontsize=7)
        axis.set_xlabel(xlabel)
        figure.tight_layout()
        png, pdf = outdir / "plots" / f"{name}.png", outdir / "plots" / f"{name}.pdf"
        figure.savefig(png, dpi=160)
        figure.savefig(pdf)
        plt.close(figure)
        figures.append({"caption": xlabel, "png": str(png.relative_to(outdir)),
                        "pdf": str(pdf.relative_to(outdir)), "image": base64.b64encode(png.read_bytes()).decode("ascii")})
    return figures


def trace_summaries(search_root, outdir):
    summaries = []
    for path in sorted(search_root.rglob("*trace.tsv")):
        if path.is_relative_to(outdir):
            continue
        try:
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if not {"name", "status"}.issubset(reader.fieldnames or []):
                    continue
                rows = list(reader)
            counts = defaultdict(int)
            for row in rows:
                counts[row["status"]] += 1
            summaries.append({"path": str(path), "recorded_tasks": len(rows), "task_statuses": dict(counts),
                              "failed_tasks": [row["name"] for row in rows if row["status"] in ("FAILED", "ABORTED")]})
        except OSError as error:
            summaries.append({"path": str(path), "error": str(error)})
    return summaries


def design_candidates(mztab, search_root, override):
    if override:
        return [[override]]
    levels = []
    folder = mztab.parent
    while True:
        candidates = sorted(set(folder.glob("*openms_design.tsv")) | set(folder.glob("*experimental_design.tsv")) |
                            set(folder.glob("ms1_design.tsv")))
        if candidates:
            levels.append(candidates)
        if folder == search_root or folder.parent == folder:
            break
        folder = folder.parent
    return levels


def find_export(mztab, consensus, candidates):
    failures = []
    for level in candidates:
        valid = []
        for design in level:
            try:
                peptides, proteins = exporter.export_ratios(mztab, consensus, design)
                valid.append((design, peptides, proteins))
            except ERRORS as error:
                failures.append(f"{design}: {error}")
        if valid:
            if len({fingerprint(item[0])["sha256"] for item in valid}) > 1:
                raise ValueError("Multiple different designs match the result; use --design")
            return valid[0]
    raise ValueError("No matching design: " + ("; ".join(failures) if failures else "design TSV absent"))


def analyze_bundle(mztab, search_root, outdir, profiles, args):
    outdir.mkdir(parents=True, exist_ok=True)
    result = {"mztab": str(mztab.resolve()), "status": "partial", "diagnostics": [],
              "identification": None, "processing": None, "comparisons": [], "benchmark": [],
              "precision_summary": [], "correlations": [], "effects": [], "answers": [], "inputs": {}}
    figures = []
    try:
        result["inputs"]["mztab"] = fingerprint(mztab)
        result["identification"], result["software"] = qc.read_identifications(mztab)
    except ERRORS as error:
        result["diagnostics"].append(f"Identification summary unavailable: {error}")
    consensus = args.consensusxml or mztab.with_suffix(".consensusXML")
    try:
        if not consensus.is_file():
            raise ValueError(f"consensusXML absent: {consensus}; executed reference and protein ratios cannot be established")
        design, peptides, proteins = find_export(mztab, consensus, design_candidates(mztab, search_root, args.design))
        result["inputs"].update({"consensusxml": fingerprint(consensus), "design": fingerprint(design)})
        groups, _ = exporter.read_design(design)
        reference, chemistry, normalized = exporter.read_settings(ET.parse(consensus).getroot())
        labels = [f"[{label}]" for label in chemistry]
        result["reference_channel"] = reference
        result["label_definitions"] = labels
        result["normalization_enabled"] = normalized
        for filename, fields, rows in (("peptide_ratios.tsv", exporter.PEPTIDE_FIELDS, peptides),
                                       ("protein_ratios.tsv", exporter.PROTEIN_FIELDS, proteins)):
            exporter.write_table(outdir / filename, fields, rows)
        result["comparisons"], grouped = qc.summarize_ratios(text_rows(peptides), text_rows(proteins))
        table(outdir / "summary.tsv", result["comparisons"], qc.SUMMARY_FIELDS)
        result["status"] = "analyzed"
        try:
            result["processing"] = qc.read_processing(consensus)
        except ERRORS as error:
            result["diagnostics"].append(f"Final feature summary unavailable: {error}")
        runs = {run for group in groups.values() for _, run in group["runs"]}
        name, profile = choose_profile(profiles, args.dataset, mztab, runs)
        result["dataset"] = name
        result["sources"] = profile.get("sources", [])
        result["questions"] = profile.get("questions", [])
        result["profile"] = profile
        unknown = sorted(runs - profile.get("runs", {}).keys())
        if unknown:
            result["diagnostics"].append("No scientific run mapping for: " + ", ".join(unknown))
        all_precision, all_correlations = [], []
        protein_pairs = []
        try:
            for level, rows in (("peptide", peptides), ("protein", proteins)):
                pairs = list(pair_ratios(rows, groups, labels, reference))
                if level == "protein":
                    protein_pairs = pairs
                result["benchmark"].extend(benchmark(pairs, groups, profile, level))
                precision, correlations = reproducibility(pairs, groups, profile, level)
                all_precision.extend(precision)
                all_correlations.extend(correlations)
            result["effects"] = biological_effects(protein_pairs, groups, profile)
        except ERRORS as error:
            result["diagnostics"].append(f"Scientific comparison unavailable: {error}")
        result["correlations"] = all_correlations
        precision_groups = defaultdict(list)
        for row in all_precision:
            key = tuple(row[field] for field in PRECISION_FIELDS[:6])
            precision_groups[key].append(row)
        for key, rows in sorted(precision_groups.items()):
            complete = [row["cv_percent"] for row in rows if row["complete"] and row["cv_percent"] is not None]
            partial = [row["cv_percent"] for row in rows if row["cv_percent"] is not None]
            result["precision_summary"].append(dict(zip(PRECISION_FIELDS[:6], key)) |
                                               {"entities": len(rows), "complete_entities": sum(row["complete"] for row in rows),
                                                "entities_with_cv": len(partial),
                                                "median_cv_percent_complete": statistics.median(complete) if complete else None,
                                                "median_cv_percent_available": statistics.median(partial) if partial else None})
        table(outdir / "benchmark.tsv", result["benchmark"], BENCHMARK_FIELDS)
        table(outdir / "replicate_precision.tsv", all_precision, PRECISION_FIELDS)
        table(outdir / "replicate_correlations.tsv", all_correlations, CORRELATION_FIELDS)
        table(outdir / "biological_effects.tsv", result["effects"], EFFECT_FIELDS)
        if not args.no_plots:
            try:
                figures = qc.make_plots(outdir, result["comparisons"], grouped, [])
                figures.extend(scientific_plots(outdir, result))
            except (ImportError, RuntimeError, ValueError, OSError) as error:
                result["diagnostics"].append(f"Plots omitted: {error}; tables and report remain available")
    except ERRORS as error:
        result["diagnostics"].append(f"Stored-ratio analysis unavailable: {error}")
    result["answers"] = scientific_answers(result)
    result["plots"] = [{key: value for key, value in figure.items() if key != "image"} for figure in figures]
    (outdir / "metrics.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result, figures


def scientific_answers(result):
    answers = []
    identification = result["identification"]
    if identification:
        answers.append(f"Identification: {identification['unique_resolved_psms']} distinct resolved PSMs and "
                       f"{identification['distinct_modified_peptides']} modified peptides in mzTab; "
                       "PSM scores and decoy annotations are recorded in metrics.json.")
    for row in result["comparisons"]:
        answers.append(f"Coverage: {row['level']}, group {row['fraction_group']}, channel "
                       f"{row['numerator_channel']}/{row['reference_channel']}, {row['normalization']}: "
                       f"{row['reported']}/{row['total_entities']} exported entities have a positive finite ratio.")
    for row in result["benchmark"]:
        answers.append(f"Accuracy: {row['level']} {row['peptide_class']}, {row['runs']}, "
                       f"{row['numerator_label']}/{row['denominator_label']}, {row['normalization']} "
                       f"({row['estimate_kind']}): expected {qc.format_value(row['expected_ratio'])}, "
                       f"median {qc.format_value(row['median_ratio'])}, bias {qc.format_value(row['median_bias_percent'])}%, "
                       f"log2 IQR {qc.format_value(row['iqr_log2'])}; n={row['reported']}.")
    if not result["benchmark"]:
        answers.append("Accuracy against ground truth is unavailable without a matching mixture definition; a ratio centered on one does not establish accuracy.")
    for row in result["precision_summary"]:
        answers.append(f"Reproducibility: {row['level']}, {row['replicate_set']} ({row['replicate_type']}), "
                       f"{row['numerator_label']}/{row['denominator_label']}, {row['normalization']}: "
                       f"median CV {qc.format_value(row['median_cv_percent_complete'])}% among "
                       f"{row['complete_entities']} entities observed in every planned replicate.")
    if not result["precision_summary"]:
        answers.append("Replicate reproducibility is unavailable without explicit replicate-set assignments.")
    if not result["effects"]:
        answers.append("Biological treatment effects are not established: configure a paired condition contrast and biological replicate identities. Isotope channels and technical injections are not biological replicates.")
    else:
        tested = sum(row["test_status"] == "tested" for row in result["effects"])
        answers.append(f"Paired biological effects: {len(result['effects'])} protein-group estimates, {tested} tested. "
                       "Inspect biological_effects.tsv for direction, biological sample counts, and BH-adjusted p-values; raw and normalized families are separate.")
    if result.get("dataset") == "PXD039580":
        answers.append("These files evaluate the HeLa mixing benchmark, not IFNγ response. A DDA-only reanalysis cannot establish a DDA-versus-DIA advantage; matched DIA/reference outputs are required.")
    return answers


def render_report(metrics, figures):
    esc = lambda value: html.escape(str(value))
    body = ["<h1>Completed SILAC run analysis</h1>",
            "<p>Stored OpenMS ratios; no imputation or additional normalization. "
            "Missing analyses and incompatible files are listed explicitly.</p>",
            '<p><a href="metrics.json">All metrics and provenance</a></p>']
    body.extend(f"<p class='note'>{esc(note)}</p>" for note in metrics["diagnostics"])
    if metrics.get("traces"):
        body.append("<h2>Execution traces</h2><p>Recorded tasks may omit running or pending jobs; traces alone do not certify completion.</p>")
        body.append(qc.html_table(metrics["traces"], [("path", "Trace"), ("recorded_tasks", "Recorded tasks"),
                                                     ("task_statuses", "Statuses"), ("failed_tasks", "Failed tasks")]))
    for bundle in metrics["bundles"]:
        folder = bundle["output_folder"]
        body.append(f"<h2>{esc(bundle['mztab'])}</h2><p>Status: {esc(bundle['status'])}; dataset: {esc(bundle.get('dataset', 'unknown'))}</p>")
        body.extend(f"<p class='note'>{esc(note)}</p>" for note in bundle["diagnostics"])
        body.extend(f"<p><strong>Study question:</strong> {esc(question)}</p>" for question in bundle.get("questions", []))
        body.append("<ul>" + "".join(f"<li>{esc(answer)}</li>" for answer in bundle["answers"]) + "</ul>")
        links = ["metrics.json"]
        if bundle["comparisons"]:
            links += ["peptide_ratios.tsv", "protein_ratios.tsv", "summary.tsv", "benchmark.tsv",
                      "replicate_precision.tsv", "replicate_correlations.tsv", "biological_effects.tsv"]
        body.append("<p>" + " · ".join(f'<a href="{esc(folder)}/{name}">{name}</a>' for name in links) + "</p>")
        if bundle["identification"]:
            body.append(qc.html_table(bundle["identification"]["per_run"], [("run", "Run"), ("unique_psms", "Distinct PSMs"), ("psm_rows", "PSM rows")]))
        for title, key, fields in (("Ratio coverage", "comparisons", qc.SUMMARY_FIELDS),
                                   ("Mixture accuracy", "benchmark", BENCHMARK_FIELDS),
                                   ("Replicate precision", "precision_summary", []),
                                   ("Replicate agreement", "correlations", CORRELATION_FIELDS)):
            rows = bundle[key]
            if rows:
                fields = fields or list(rows[0])
                body.append(f"<h3>{title}</h3>" + qc.html_table(rows, [(field, field.replace('_', ' ')) for field in fields]))
        for source in bundle.get("sources", []):
            if source.startswith(("https://", "http://")):
                body.append(f'<p>Source: <a href="{esc(source)}">{esc(source)}</a></p>')
    for figure in figures:
        body.append(f'<figure><img alt="{esc(figure["caption"])}" src="data:image/png;base64,{figure["image"]}"><figcaption>{esc(figure["caption"])}</figcaption></figure>')
    body.append("<h2>Interpretation and methods</h2><ul>" + "".join(f"<li>{esc(note)}</li>" for note in metrics["methods"]) + "</ul>")
    return '''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>SILAC run analysis</title>
<style>body{font:15px/1.6 system-ui;background:#f5f7fa;color:#243441;margin:30px auto;padding:0 24px;max-width:1400px}
h1,h2,h3{color:#173a4c;overflow-wrap:anywhere}.note{background:#fff0d5;padding:12px}.table-scroll{overflow:auto}
table{border-collapse:collapse;background:white}td,th{padding:8px;border:1px solid #dbe4e8;white-space:nowrap}
img{max-width:100%}a{color:#176a83}li{margin-bottom:8px}</style></head><body>''' + "\n".join(body) + "</body></html>"


METHODS = [
    "Each mzTab/consensusXML/design bundle is analyzed independently. Separate cohorts or pipeline invocations are never pooled as replicates.",
    "Coverage denominator: all experiment-wide exported entities in each fraction-group/channel comparison; missing estimates are not evidence of biological absence.",
    "Raw and OpenMS-normalized ratios are separate. Unequal mixture benchmarks should primarily be assessed using raw ratios.",
    "Protein identities retain full inferred-group membership. Stored peptide evidence counts are not biological replicate counts.",
    "All channel pairs use stored ratios. Derived non-reference ratios are quotients of aggregated estimates; they can differ from a native OpenMS run with the new reference, particularly for proteins. They require both source ratios.",
    "Replicate CV is sample SD / mean on linear ratios; SD and correlations use log2 ratios. Complete-case and available-case CV summaries are separate. Technical precision does not measure biological variability.",
    "Correlation uses shared quantified entities, requires at least three, and is undefined for constant vectors. Correlation alone does not establish agreement or accuracy.",
    "Biological tests require an explicit paired condition map. Technical repeats collapse by median log2 ratio within biological replicate. At least three biological replicates and nonzero variance are required for a two-sided one-sample t-test against zero; BH adjustment is per contrast and normalization within a bundle.",
    "The paired t-test assumes independent biological replicates and approximately normal log2 effects. It is exploratory, without empirical Bayes moderation, covariates, pathway enrichment, or site localization. It does not reproduce the publications' full statistical pipelines.",
    "Input files have SHA-256 checksums. Scores are called q-values only when the mzTab metadata identifies them as such; reported q-values do not independently measure true FDR.",
    "A partial or absent result produces a diagnostic report, not a scientific conclusion or a successful-run certification. No universal scientific pass/fail threshold is imposed.",
]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="Published result directory (searched recursively), or one mzTab")
    parser.add_argument("--outdir", type=Path, help="Default: RESULTS/analysis (or beside the mzTab)")
    parser.add_argument("--design", type=Path, help="Explicit design; allowed only for one mzTab")
    parser.add_argument("--consensusxml", type=Path, help="Explicit consensusXML; allowed only for one mzTab")
    parser.add_argument("--dataset", default="auto", help="auto, none, or a dataset key in the profile JSON")
    parser.add_argument("--profiles", type=Path, default=ROOT / "conf/analysis_profiles.json")
    parser.add_argument("--no-plots", action="store_true", help="Use only standard-library report generation")
    parser.add_argument("--strict", action="store_true", help="Exit 2 after writing report if any analysis has diagnostics")
    args = parser.parse_args(argv)
    try:
        if not args.results.exists():
            raise ValueError(f"Input does not exist: {args.results}")
        search_root = args.results.resolve() if args.results.is_dir() else args.results.resolve().parent
        outdir = (args.outdir or search_root / "analysis").resolve()
        if args.results.is_file():
            mztabs = [args.results.resolve()]
        else:
            mztabs = sorted(path for path in search_root.rglob("*") if path.is_file() and path.suffix.lower() == ".mztab"
                            and not path.is_relative_to(outdir))
        if len(mztabs) != 1 and (args.design or args.consensusxml):
            raise ValueError("Explicit design/consensusXML requires exactly one mzTab; analyze that file directly")
        profiles = load_profiles(args.profiles)
        if args.dataset not in ("auto", "none") and args.dataset not in profiles:
            raise ValueError(f"Unknown dataset: {args.dataset}")
        outdir.mkdir(parents=True, exist_ok=True)
        metrics = {"schema_version": 1, "script": fingerprint(Path(__file__)), "python_version": sys.version,
                   "profiles_input": fingerprint(args.profiles), "bundles": [], "diagnostics": [], "methods": METHODS}
        metrics["analysis_helpers"] = {name: fingerprint(ROOT / "bin" / name)
                                       for name in ("export_silac_ratios.py", "silac_qc.py")}
        metrics["optional_packages"] = {}
        for package in ("matplotlib", "scipy"):
            try:
                metrics["optional_packages"][package] = importlib.metadata.version(package)
            except importlib.metadata.PackageNotFoundError:
                metrics["optional_packages"][package] = None
        metrics["traces"] = trace_summaries(search_root, outdir)
        if args.results.is_dir():
            for design in sorted(search_root.rglob("*openms_design.tsv")):
                if not design.is_relative_to(outdir) and not any(design.parent in path.parents for path in mztabs):
                    metrics["diagnostics"].append(f"No mzTab output beneath metadata directory: {design.parent}. Quantitative analysis is unavailable for this run.")
        figures = []
        if not mztabs:
            metrics["diagnostics"].append("No mzTab outputs found. The run may be unfinished or outputs unpublished; quantitative analysis is unavailable.")
        for index, mztab in enumerate(mztabs, 1):
            folder = f"bundle_{index:03d}_" + hashlib.sha256(str(mztab).encode()).hexdigest()[:8]
            result, plots = analyze_bundle(mztab, search_root, outdir / folder, profiles, args)
            result["output_folder"] = folder
            metrics["bundles"].append(result)
            figures.extend(plots)
            print(f"{mztab}: {result['status']} ({len(result['diagnostics'])} diagnostics)")
        (outdir / "metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n")
        (outdir / "report.html").write_text(render_report(metrics, figures))
        print(f"Analysis report: {outdir / 'report.html'}")
        return 2 if args.strict and (metrics["diagnostics"] or any(row["diagnostics"] for row in metrics["bundles"])) else 0
    except ERRORS as error:
        parser.exit(1, f"Cannot analyze completed run: {error}\n")


if __name__ == "__main__":
    sys.exit(main())
