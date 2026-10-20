#!/usr/bin/env python3
"""Compare retained SILAC ratios with deposited measurements and known mixtures.

Requires numpy and matplotlib (available in the local silac conda environment).
Reads existing tables only; no search, quantification, or normalization is run.
"""

import argparse
from collections import defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import textwrap
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
H1975 = ROOT / "results/pxd001101_h1975_br1_pilot_qpx_20261009/PXD001101_H1975_BR1_pilot/analysis/bundle_001_0ae98b1f/peptide_ratios.tsv"
ARCHIVE = ROOT.parent / "sdrf_annotations/datasets/PXD001101/input/H1975_MaxQuantresults.zip"
HELA = ROOT / "results/pxd039580_raw_all6_qpx_20261007_131122/analysis"
ECOLI = ROOT / "results/PXD003327_duplex/tables/peptide_ratios.tsv"
TRIPLEX = ROOT / "results/analysis_repeat_20261014_163325/ecoli_triplex_historical/bundle_001_94134937/peptide_ratios.tsv"
BR = "H1975_pST_BR1"
BLUE, ORANGE, GREEN, GRAY = "#247ba0", "#d17b32", "#39856a", "#657085"
SOURCES = {
    "PXD001101_paper": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6461560/",
    "PXD001101_results": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/01/PXD001101/H1975_MaxQuantresults.zip",
    "PXD039580_paper": "https://www.nature.com/articles/s41467-023-43919-3",
    "PXD039580_supplement": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-43919-3/MediaObjects/41467_2023_43919_MOESM6_ESM.xlsx",
    "Ecoli_mixtures_Table_I": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4974355/",
}


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(out, name, rows):
    assert rows, name
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (out / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def positive(value):
    try:
        number = float(value)
        return number if math.isfinite(number) and number > 0 else None
    except (ValueError, TypeError):
        return None


def clean(row):
    return row.get("Reverse") != "+" and row.get("Contaminant") != "+"


def fingerprint(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return {"path": str(path.resolve()), "sha256": h.hexdigest()}


def canonical(form):
    return form.strip("_").replace("(ph)", "(Phospho)").replace("(ox)", "(Oxidation)").replace("(Carbamidomethyl)", "")


def bare(form):
    return re.sub(r"\([^)]*\)", "", form)


def quantitative(rows):
    """One native ratio per run-group, entity, channel; retain only raw values."""
    values = defaultdict(dict)
    for row in rows:
        value = positive(row["ratio"])
        if row["normalization"] != "raw" or value is None:
            continue
        key = row["runs"], row.get("peptidoform", row.get("protein_group"))
        channel = int(row["numerator_channel"])
        assert channel not in values[key], (key, channel)
        values[key][channel] = value
    return values


def format_ratio(x, _=None):
    return f"{2.0 ** x:.3g}"


def ratio_axis(axis, direction="both"):
    if direction in ("both", "x"):
        axis.xaxis.set_major_formatter(FuncFormatter(format_ratio))
    if direction in ("both", "y"):
        axis.yaxis.set_major_formatter(FuncFormatter(format_ratio))


def finish(fig, title, subtitle, note):
    fig.suptitle(title, x=0.08, ha="left", fontsize=17, fontweight="bold", y=0.985)
    fig.text(0.08, 0.935, subtitle, fontsize=10, color=GRAY)
    note = "\n".join(textwrap.fill(line, width=int(fig.get_figwidth() * 13.2)) for line in note.splitlines())
    fig.text(0.08, 0.025, note, fontsize=9, color=GRAY, va="bottom")
    fig.tight_layout(rect=(0.025, 0.115, 0.99, 0.90), h_pad=3, w_pad=3)


def scatter(axis, rows, title):
    assert len(rows) >= 3, title
    xx = np.log2([r["reference_ratio"] for r in rows])
    yy = np.log2([r["our_ratio"] for r in rows])
    lo, hi = min(xx.min(), yy.min()), max(xx.max(), yy.max())
    pad = max(0.12, (hi - lo) * .06)
    lo, hi = lo - pad, hi + pad
    axis.scatter(xx, yy, s=10, alpha=.4, color=BLUE, edgecolors="none", rasterized=True)
    axis.plot([lo, hi], [lo, hi], color=GRAY, lw=1, ls="--")
    axis.set(xlim=(lo, hi), ylim=(lo, hi), title=title,
             xlabel="Deposited ratio (log2-spaced)", ylabel="Our raw ratio (log2-spaced)")
    ratio_axis(axis)
    r = float(np.corrcoef(xx, yy)[0, 1])
    error = float(np.median(np.abs(yy - xx)))
    axis.text(.04, .96, f"n = {len(rows):,}   Pearson r = {r:.3f}\nmedian |log2 difference| = {error:.3f}",
              transform=axis.transAxes, va="top", fontsize=9,
              bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
    return {"n": len(rows), "pearson_log2": r, "median_absolute_log2_difference": error,
            "median_log2_difference": float(np.median(yy - xx))}


def boxes(axis, arrays, labels, colors, expected=None):
    logs = [np.log2(a) for a in arrays]
    result = axis.boxplot(logs, tick_labels=labels, patch_artist=True, widths=.55,
                          showfliers=False, whis=(5, 95), medianprops={"color": "#152b3b", "lw": 1.5})
    for patch, color in zip(result["boxes"], colors):
        patch.set(facecolor=color, alpha=.65)
    if expected is not None:
        axis.axhline(math.log2(expected), color=GREEN, ls="--", label=f"Expected = {expected:.4g}")
        axis.legend(loc="upper right", fontsize=9)
    axis.set_ylabel("Ratio (log2-spaced)")
    ratio_axis(axis, "y")
    axis.margins(y=.22)
    for i, a in enumerate(arrays, 1):
        axis.text(i, .02, f"n={len(a):,}\nmed={np.median(a):.3g}",
                  transform=axis.get_xaxis_transform(), ha="center", fontsize=8)


def h1975_figures(out, save):
    ours_rows = read(H1975)
    ours = quantitative(ours_rows)
    with zipfile.ZipFile(ARCHIVE) as archive:
        sites = [r for r in csv.DictReader(io.TextIOWrapper(archive.open("PhosphoSTYSites.txt")), delimiter="\t") if clean(r)]
        peptides = [r for r in csv.DictReader(io.TextIOWrapper(archive.open("peptides.txt")), delimiter="\t") if clean(r) and r["Experiment " + BR]]
    sites = [r for r in sites if (positive(r["Localization prob " + BR]) or 0) > .75]
    forms_by_seq = defaultdict(set)
    for _, form in ours:
        forms_by_seq[bare(canonical(form))].add(canonical(form))
    # All BR1 class-I site rows count toward ambiguity, even if a ratio is absent.
    sites_by_seq = defaultdict(list)
    for site in sites:
        sites_by_seq[bare(canonical(site["Modified sequence"]))].append(site)
    reference_peptides = {r["Sequence"]: r for r in peptides}
    assert len(reference_peptides) == len(peptides)
    matched, distributions = [], []
    for contrast, channel in [("M/L", 1), ("H/M", 3)]:
        for site in sites:
            for normalization, suffix in [("raw", ""), ("MaxQuant normalized", " normalized")]:
                value = positive(site["Ratio " + contrast + suffix + " " + BR])
                if value:
                    distributions.append(dict(population="Deposited BR1 class-I sites", normalization=normalization,
                                              contrast=contrast, entity=site["id"], ratio=value))
        for (_, form), channels in sorted(ours.items()):
            if channel not in channels:
                continue
            value = 1 / channels[channel] if channel == 1 else channels[channel]
            normalized_form = canonical(form)
            seq = bare(normalized_form)
            if "(Phospho)" in form:
                distributions.append(dict(population="Our BR1 phosphopeptide forms", normalization="raw",
                                          contrast=contrast, entity=form, ratio=value))
            reference, population = None, None
            only_one_form = len(forms_by_seq[seq]) == 1
            pep = reference_peptides.get(seq)
            if only_one_form and "(" not in normalized_form and pep and not pep["Phospho (STY) site IDs"] and not pep["Oxidation (M) site IDs"]:
                reference, population = pep, "No variable PTM"
            candidates = sites_by_seq[seq]
            if only_one_form and len(candidates) == 1 and normalized_form.count("(Phospho)") == 1:
                candidate = candidates[0]
                if normalized_form == canonical(candidate["Modified sequence"]) and "(" not in normalized_form.replace("(Phospho)", ""):
                    reference, population = candidate, "Assigned monophosphopeptide"
            if reference is not None:
                reference_ratio = positive(reference["Ratio " + contrast + " " + BR])
                if reference_ratio:
                    matched.append(dict(population=population, sequence=seq, our_peptidoform=form,
                                        contrast=contrast, our_ratio=value, reference_ratio=reference_ratio,
                                        reference_normalized_ratio=reference["Ratio " + contrast + " normalized " + BR],
                                        reference_id=reference["id"], our_feature_count=next(r["ratio_count"] for r in ours_rows if r["peptidoform"] == form and int(r["numerator_channel"]) == channel),
                                        log2_difference=math.log2(value / reference_ratio)))
    write(out, "pxd001101_distributions.tsv", distributions)
    write(out, "pxd001101_matched_ratios.tsv", matched)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8))
    summary = []
    populations = [("Our BR1 phosphopeptide forms", "raw", BLUE, "Our raw forms"),
                   ("Deposited BR1 class-I sites", "raw", ORANGE, "Deposited BR1 raw sites"),
                   ("Deposited BR1 class-I sites", "MaxQuant normalized", GREEN, "Deposited BR1 normalized sites")]
    for axis, contrast in zip(axes, ["M/L", "H/M"]):
        selected = [[r["ratio"] for r in distributions if r["population"] == p and r["normalization"] == n and r["contrast"] == contrast]
                    for p, n, _, _ in populations]
        # ECDF includes the entire distribution, with no histogram truncation.
        for values, (p, n, color, label) in zip(selected, populations):
            x = np.sort(np.log2(values))
            axis.step(x, np.arange(1, len(x) + 1) / len(x), color=color,
                      label=f"{label}\nn={len(x):,}; median={np.median(values):.3f}", lw=1.8)
            summary.append(dict(dataset="PXD001101", contrast=contrast, population=p, normalization=n,
                                n=len(values), median=float(np.median(values))))
        axis.axvline(0, color=GRAY, lw=1, ls="--")
        axis.set(xlabel=f"{contrast} (log2-spaced)", ylabel="Cumulative fraction", title=contrast)
        ratio_axis(axis, "x")
        axis.legend(fontsize=8, loc="lower right")
    finish(fig, "PXD001101 · EGF and erlotinib responses", "H1975 BR1: M/L = EGF/control; H/M = erlotinib+EGF/EGF",
           "Descriptive populations: unlocalized peptide forms versus localized site rows; these are different units.\n"
           "One biological replicate. Dashed ratio 1 means unchanged, not a universal expected phosphopeptide ratio.")
    save(fig, "pxd001101_distributions")
    metrics = {}
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for i, population in enumerate(["Assigned monophosphopeptide", "No variable PTM"]):
        for j, contrast in enumerate(["M/L", "H/M"]):
            rr = [r for r in matched if r["population"] == population and r["contrast"] == contrast]
            metrics[population + " " + contrast] = scatter(axes[i, j], rr, population + " · " + contrast)
    finish(fig, "PXD001101 · Ratios on matched peptide candidates", "Deposited BR1 raw versus our QPX-enabled BR1 raw output; dashed line = exact agreement",
           "Top: identical assigned mono-phospho form, one our form/sequence, one deposited class-I site row/sequence.\n"
           "Our localization remains unvalidated; charge/evidence may differ. Bottom: one form/sequence, no variable PTMs or reference site links.")
    save(fig, "pxd001101_matched_ratios")
    # Preselected signaling examples from the earlier study review, with explicit missingness.
    targets = [("MAPK1", "T185", "VADPDHDHTGFLT(Phospho)EY(Phospho)VATR"),
               ("EGFR", "T693", "ELVEPLT(Phospho)PSGEAPNQ(Deamidated)ALLR"),
               ("MAPK3", "T202", "IADPEHDHTGFLT(Phospho)EY(Phospho)VATR")]
    signaling = []
    for gene, site_name, form in targets:
        ref = [r for r in sites if r["Amino acid"] + r["Position"] == site_name and bare(canonical(r["Modified sequence"])) == bare(form)]
        assert len(ref) == 1, (gene, site_name)
        ref = ref[0]
        channels = next((v for (_, f), v in ours.items() if f == form), {})
        for contrast, channel in [("M/L", 1), ("H/M", 3)]:
            signaling.append(dict(gene=gene, site_candidate=site_name, contrast=contrast,
                                  reference_BR1_raw=float(ref["Ratio " + contrast + " " + BR]),
                                  reference_BR1_normalized=float(ref["Ratio " + contrast + " normalized " + BR]),
                                  reference_pooled_normalized=float(ref["Ratio " + contrast + " normalized"]),
                                  our_raw=(1 / channels[channel] if channel == 1 else channels[channel]) if channel in channels else None,
                                  our_form=form, reference_best_form=ref["Modified sequence"],
                                  identity_status="sequence-linked candidate; site/form equivalence unvalidated"))
    write(out, "pxd001101_signaling_candidates.tsv", signaling)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for axis, contrast in zip(axes, ["M/L", "H/M"]):
        rr = [r for r in signaling if r["contrast"] == contrast]
        for offset, key, color, label in [(-.27, "reference_pooled_normalized", GRAY, "Paper pooled normalized"),
                                           (-.09, "reference_BR1_normalized", GREEN, "Deposited BR1 normalized"),
                                           (.09, "reference_BR1_raw", ORANGE, "Deposited BR1 raw"),
                                           (.27, "our_raw", BLUE, "Our BR1 raw")]:
            x = np.array([i + offset for i, r in enumerate(rr) if positive(r[key])])
            y = np.array([r[key] for r in rr if positive(r[key])])
            axis.scatter(x, np.log2(y), color=color, s=55, label=label, zorder=3)
            for xx, yy in zip(x, y):
                axis.annotate(f"{yy:.2f}", (xx, math.log2(yy)), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=7)
        axis.text(2.27, -.15, "No ratio", rotation=90, ha="center", fontsize=8, color=BLUE)
        axis.axhline(0, color=GRAY, ls="--", lw=1)
        axis.set(xticks=range(3), xticklabels=["MAPK1\nT185-linked", "EGFR\nT693-linked", "MAPK3\nT202-linked"],
                 ylabel=f"{contrast} (log2-spaced)", title=contrast, xlim=(-.55, 2.55))
        ratio_axis(axis, "y")
        axis.margins(y=.30)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, loc="upper center", bbox_to_anchor=(.53, .89), ncol=4, frameon=False)
    finish(fig, "PXD001101 · Selected signaling examples", "Replicate and normalization explain why the paper's reported values differ from raw BR1",
           "Sequence-linked examples, not validated site matches. MAPK1: doubly phosphorylated form; EGFR: our form also has deamidation.\n"
           "MAPK3 has no retained quantity. Paper pooled values summarize multiple biological experiments; our pilot covers BR1 only.")
    save(fig, "pxd001101_signaling_examples")
    return summary, metrics


def accession_set(value):
    return ";".join(sorted({p.split("|")[1] if "|" in p else p for p in value.split(";")}))


def hela_figures(out, save):
    ours = quantitative(read(HELA / "bundle_001_e87c1a55/protein_ratios.tsv"))
    measurements = []
    for (run, group), channels in ours.items():
        assert "A3R1_" in run or "A3R2_" in run
        mixture = "Mix1" if "A3R1_" in run else "Mix2"
        for contrast, value in [("M/L", channels.get(2)), ("H/M", channels[3] / channels[2] if 2 in channels and 3 in channels else None)]:
            if value is not None:
                measurements.append(dict(method="Ours", mixture=mixture, run=run, protein_group=accession_set(group),
                                         contrast=contrast, ratio=value))
    for row in read(HELA / "publication_dda_ratios.tsv"):
        if row["ratio"] not in ("M/L", "H/M") or not positive(row["num_ratio"]):
            continue
        assert re.fullmatch(r"AAAR[12]_R[123]", row["Experiment"])
        measurements.append(dict(method="Published DDA", mixture="Mix" + row["Experiment"][4],
                                 run=row["Experiment"], protein_group=accession_set(row["Protein.Ids"]),
                                 contrast=row["ratio"], ratio=float(row["num_ratio"])))
    grouped = defaultdict(list)
    seen = set()
    for row in measurements:
        identity = row["method"], row["run"], row["protein_group"], row["contrast"]
        assert identity not in seen, identity
        seen.add(identity)
        grouped[(row["method"], row["mixture"], row["contrast"], row["protein_group"])].append(row["ratio"])
    means = {key: float(np.mean(v)) for key, v in grouped.items()}
    write(out, "pxd039580_per_injection_ratios.tsv", measurements)
    write(out, "pxd039580_protein_means.tsv", [dict(method=k[0], mixture=k[1], contrast=k[2], protein_group=k[3],
                                                  ratio_mean=v, injections=len(grouped[k])) for k, v in sorted(means.items())])
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    paired, summaries = [], []
    for i, mixture in enumerate(["Mix1", "Mix2"]):
        for j, contrast in enumerate(["M/L", "H/M"]):
            reference = {k[3]: v for k, v in means.items() if k[:3] == ("Published DDA", mixture, contrast)}
            ours = {k[3]: v for k, v in means.items() if k[:3] == ("Ours", mixture, contrast)}
            shared = sorted(reference.keys() & ours.keys())
            expected = 15 / 70 if mixture == "Mix1" and contrast == "M/L" else 2 if contrast == "M/L" else 1
            arrays = [list(reference.values()), [reference[k] for k in shared], [ours[k] for k in shared], list(ours.values())]
            boxes(axes[i, j], arrays, ["Paper\nall", "Paper\nshared", "Ours\nshared", "Ours\nall"],
                  [ORANGE, ORANGE, BLUE, BLUE], expected)
            axes[i, j].set_title(mixture + " · " + contrast)
            for name, values in zip(["Paper all", "Paper shared", "Ours shared", "Ours all"], arrays):
                summaries.append(dict(dataset="PXD039580", contrast=contrast, population=mixture + " " + name,
                                      normalization="published column provenance unresolved" if name.startswith("Paper") else "raw",
                                      n=len(values), median=float(np.median(values)), expected=expected,
                                      median_bias_percent=100 * (float(np.median(values)) / expected - 1)))
            paired.extend(dict(mixture=mixture, contrast=contrast, protein_group=k, reference_ratio=reference[k], our_ratio=ours[k]) for k in shared)
    finish(fig, "PXD039580 · HeLa expected mixtures versus measured ratios", "Protein-group linear mean across available injections (1–3); each group contributes once per panel",
           "Boxes: median and middle 50%; whiskers: 5th–95th percentile; outliers hidden but included in statistics.\n"
           "Shared = identical accession sets. Our H/M = stored H/L ÷ M/L per injection; published H/M is native. Published normalization is unresolved.")
    save(fig, "pxd039580_expected_ratios")
    write(out, "pxd039580_matched_ratios.tsv", paired)
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    metrics = {}
    for i, mixture in enumerate(["Mix1", "Mix2"]):
        for j, contrast in enumerate(["M/L", "H/M"]):
            rr = [r for r in paired if r["mixture"] == mixture and r["contrast"] == contrast]
            metrics[mixture + " " + contrast] = scatter(axes[i, j], rr, mixture + " · " + contrast)
    finish(fig, "PXD039580 · Published DDA versus our protein ratios", "Matched accession sets; one available-injection linear mean per group; dashed line = exact agreement",
           "Three technical injections per mixture, not biological replicates. Shared groups can have different peptide support or injection availability.\n"
           "Our H/M is derived per injection. Reference normalization provenance is unresolved; agreement does not by itself establish mixture accuracy.")
    save(fig, "pxd039580_matched_ratios")
    return summaries, metrics


def ecoli_figures(out, save):
    observations, summaries = [], []
    for dataset, path in [("Duplex Nextflow", ECOLI), ("Triplex historical Bash", TRIPLEX)]:
        rows = read(path)
        lookup = {(r["runs"], r["peptidoform"]): r["sequence"] for r in rows}
        for (run, form), channels in quantitative(rows).items():
            seq = lookup[run, form]
            kind = "Single K" if seq.count("K") == 1 and seq.count("R") == 0 else "Single R" if seq.count("R") == 1 and seq.count("K") == 0 else None
            if kind is None:
                continue
            if dataset == "Duplex Nextflow":
                expected = 1 if "1-1" in run else .25 if kind == "Single K" else 4
                contrasts = [("H/L", channels.get(2), expected)]
            else:
                contrasts = [("M/L", channels.get(2), 2 if kind == "Single K" else .5),
                             ("H/M", channels[3] / channels[2] if 2 in channels and 3 in channels else None,
                              2 if kind == "Single K" else .5)]
            for contrast, value, expected in contrasts:
                if value is not None:
                    observations.append(dict(dataset=dataset, run=run, peptide_class=kind, peptidoform=form,
                                             contrast=contrast, ratio=value, expected=expected,
                                             estimate_kind="derived H/L divided by M/L" if contrast == "H/M" else "native"))
    write(out, "ecoli_ratios.tsv", observations)
    for dataset, name in [("Duplex Nextflow", "pxd003327_duplex_expected_ratios"), ("Triplex historical Bash", "pxd003327_historical_triplex_expected_ratios")]:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        panels = ["Chris_Ecoli_4-1", "Chris_Ecoli_1-1"] if dataset == "Duplex Nextflow" else ["M/L", "H/M"]
        for axis, panel in zip(axes, panels):
            arrays, labels, expected_values = [], [], []
            for kind in ["Single K", "Single R"]:
                rr = [r for r in observations if r["dataset"] == dataset and r["peptide_class"] == kind
                      and (r["run"] == panel if dataset == "Duplex Nextflow" else r["contrast"] == panel)]
                values = [r["ratio"] for r in rr]
                expected = rr[0]["expected"]
                arrays.append(values)
                expected_values.append(expected)
                labels.append(kind)
                summaries.append(dict(dataset="E. coli " + dataset, contrast=rr[0]["contrast"], population=panel + " " + kind,
                                      normalization="raw", n=len(values), median=float(np.median(values)), expected=expected,
                                      median_bias_percent=100 * (float(np.median(values)) / expected - 1)))
            boxes(axis, arrays, labels, [BLUE, ORANGE])
            for pos, expected in enumerate(expected_values, 1):
                axis.plot([pos - .35, pos + .35], [math.log2(expected)] * 2, ls="--", color=GREEN, lw=2)
                axis.annotate(f"Expected {expected:g}", (pos + .36, math.log2(expected)), fontsize=8, color=GREEN)
            axis.set_title(panel + (" · H/L" if dataset == "Duplex Nextflow" else ""))
            axis.set_xlim(.4, 2.95)
        finish(fig, "E. coli · " + dataset + " ratio recovery", "Known residue-specific mixtures from PyQuant Table I; dashed green marks = intended ratio",
               "Exactly one K/no R or one R/no K; stored raw peptide ratios. Boxes: IQR; whiskers: 5th–95th percentile.\n" +
               ("Duplex has light and heavy channels only, so H/L is the relevant comparison." if dataset == "Duplex Nextflow" else
                "Historical Bash predecessor, not a current Nextflow triplex validation. H/M is derived for each peptide as H/L ÷ M/L.") +
               " Local accession label PXD003327 has a documented PXD003845 provenance ambiguity.")
        save(fig, name)
    return summaries


def h1975_boxplot(out, save):
    rows = read(out / "pxd001101_distributions.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(9, 5.5))
    for axis, contrast in zip(axes, ["M/L", "H/M"]):
        arrays = [[float(r["ratio"]) for r in rows
                   if r["contrast"] == contrast and r["normalization"] == "raw"
                   and r["population"] == population]
                  for population in ["Deposited BR1 class-I sites", "Our BR1 phosphopeptide forms"]]
        boxes(axis, arrays, ["Paper BR1", "Ours BR1"], [ORANGE, BLUE])
        axis.axhline(0, color=GRAY, ls="--", lw=1)
        axis.set_title(contrast)
    finish(fig, "PXD001101 · Paper versus ours", "H1975 BR1 · raw SILAC ratios",
           "Paper: 2,011 localized phosphosite rows; ours: 1,536 phosphopeptide forms (localization unvalidated).\n"
           "Boxes: middle 50%; line: median; whiskers: 5th–95th percentiles. Dashed line: ratio 1.")
    save(fig, "pxd001101_boxplot")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=ROOT / "results/paper_ratio_comparisons_20261015")
    args = parser.parse_args()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 170})
    figures = []
    with PdfPages(out / "paper_vs_ours_all_datasets.pdf") as pdf:
        def save(fig, name):
            fig.savefig(out / (name + ".pdf"))
            fig.savefig(out / (name + ".png"))
            pdf.savefig(fig)
            plt.close(fig)
            figures.append(name)
        h_summary, h_metrics = h1975_figures(out, save)
        hela_summary, hela_metrics = hela_figures(out, save)
        e_summary = ecoli_figures(out, save)
        h1975_boxplot(out, save)
    write(out, "ratio_summary.tsv", h_summary + hela_summary + e_summary)
    (out / "agreement_metrics.json").write_text(json.dumps({"PXD001101": h_metrics, "PXD039580": hela_metrics}, indent=2) + "\n")
    inputs = [H1975, ARCHIVE, HELA / "bundle_001_e87c1a55/protein_ratios.tsv",
              HELA / "publication_dda_ratios.tsv", HELA / "publication_reference/Supplementary_Data_3.xlsx", ECOLI, TRIPLEX]
    provenance = {"sources": SOURCES, "inputs": [fingerprint(p) for p in inputs],
                  "script": fingerprint(Path(__file__)), "numpy": np.__version__, "matplotlib": matplotlib.__version__,
                  "figures": figures, "output_tables": [fingerprint(p) for p in sorted(out.glob("*.tsv"))]}
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"outdir": str(out), "figures": figures, "agreement": {"PXD001101": h_metrics, "PXD039580": hela_metrics}}, indent=2))


if __name__ == "__main__":
    main()
