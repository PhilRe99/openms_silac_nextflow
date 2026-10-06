#!/usr/bin/env python3
"""Export stored OpenMS peptide/protein ratios without requantifying abundances."""

import argparse
import csv
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path


COMMON_FIELDS = [
    "fraction_group", "runs", "numerator_channel", "reference_channel",
    "numerator_label_definition", "reference_label_definition",
    "numerator_sample", "reference_sample", "numerator_condition",
    "reference_condition", "normalization", "ratio", "log2_ratio",
    "ratio_count", "ratio_status",
]
PEPTIDE_FIELDS = ["sequence", "modifications", "peptidoform"] + COMMON_FIELDS
PROTEIN_FIELDS = ["protein_group"] + COMMON_FIELDS
PEPTIDE_PREFIX = "opt_global_MS1Label:peptide_ratio"


def user_params(node):
    return {child.get("name"): child.get("value") for child in node.findall("UserParam")}


def read_design(path):
    sections = [[], []]
    section = 0
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                if sections[0]:
                    section = 1
            else:
                sections[section].append(line)
    samples = {row["Sample"]: row for row in csv.DictReader(sections[1], delimiter="\t")}
    groups = {}
    for row in csv.DictReader(sections[0], delimiter="\t"):
        group_id = int(row["Fraction_Group"])
        channel = int(row["Label"])
        group = groups.setdefault(group_id, {"runs": set(), "samples": {}})
        group["runs"].add((int(row["Fraction"]), Path(row["Spectra_Filepath"]).stem))
        sample = row["Sample"]
        if sample not in samples:
            raise ValueError(f"Sample {sample} is absent from the design sample table")
        previous = group["samples"].setdefault(channel, sample)
        if previous != sample:
            raise ValueError(f"Fraction group {group_id}, channel {channel} maps to multiple samples")
    if not groups:
        raise ValueError("Experimental design contains no runs")
    return groups, samples


def read_settings(root):
    processing = [node for node in root.findall("dataProcessing")
                  if node.find("software") is not None
                  and node.find("software").get("name") == "MS1LabeledWorkflow"]
    if len(processing) != 1:
        raise ValueError("Expected one MS1LabeledWorkflow dataProcessing record")
    params = user_params(processing[0])
    reference = int(params["parameter: ratios:reference_channel"])
    labels = re.findall(r"\[([^\[\]]*)\]", params["parameter: labels"])
    if not labels or not 1 <= reference <= len(labels):
        raise ValueError("Invalid executed labels/reference channel in consensusXML")
    normalize = params["parameter: ratios:normalize"]
    if normalize not in ("true", "false"):
        raise ValueError(f"Invalid ratio normalization setting: {normalize}")
    return reference, labels, normalize == "true"


def check_design_maps(root, groups):
    expected = {(group_id, fraction, run, channel)
                for group_id, group in groups.items()
                for fraction, run in group["runs"] for channel in group["samples"]}
    actual = set()
    for column in root.findall("./mapList/map"):
        params = user_params(column)
        actual.add((int(params["fraction_group"]), int(params["fraction"]),
                    Path(column.get("name")).stem, int(params["channel_id"]) + 1))
    if actual != expected:
        raise ValueError("Design run/fraction/channel mapping does not match consensusXML")


def array(value, converter):
    if value in (None, "", "null"):
        return []
    if not value.startswith("[") or not value.endswith("]"):
        raise ValueError(f"Expected an OpenMS numeric array, got {value!r}")
    return [None if item.strip() == "null" else converter(item.strip())
            for item in value[1:-1].split(",") if item.strip()]


def read_ratio_arrays(params, prefix, channel_suffix, normalized):
    """Normalized quantities share the raw quantity's group/channel/count arrays."""
    result = []
    for normalization in (["raw", "normalized"] if normalized else ["raw"]):
        key = prefix + ("_normalized" if normalization == "normalized" else "")
        if params.get(key) in (None, "", "null"):
            continue
        values = array(params[key], float)
        groups = array(params.get(prefix + "_fraction_group"), int)
        channels = array(params.get(prefix + channel_suffix), int)
        counts = array(params.get(prefix + "_count"), int)
        if len({len(values), len(groups), len(channels), len(counts)}) != 1:
            raise ValueError(f"Ratio arrays have different lengths: {key}")
        for group, channel, value, count in zip(groups, channels, values, counts):
            if group is None or channel is None or (count is not None and count < 0):
                raise ValueError(f"Invalid ratio identity/count: {key}")
            result.append(((group, channel, normalization), (value, count)))
    return result


def add_ratio(ratios, key, observation):
    if key in ratios:
        old_value, old_count = ratios[key]
        value, count = observation
        same_value = old_value == value or (
            old_value is not None and value is not None
            and ((math.isnan(old_value) and math.isnan(value))
                 or math.isclose(old_value, value, rel_tol=1e-12, abs_tol=0.0)))
        if not same_value or old_count != count:
            raise ValueError(f"Conflicting repeated ratio: {key}")
    else:
        ratios[key] = observation


def read_peptides(path, normalized):
    peptides, ratios = {}, {}
    header = None
    with path.open(newline="") as handle:
        for line_number, record in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not record:
                continue
            if record[0] == "PEH":
                header = record
                required = {"sequence", "modifications", PEPTIDE_PREFIX,
                            PEPTIDE_PREFIX + "_channel", PEPTIDE_PREFIX + "_count",
                            PEPTIDE_PREFIX + "_fraction_group"}
                if not required.issubset(header):
                    raise ValueError(f"Missing mzTab peptide fields: {sorted(required - set(header))}")
            elif record[0] == "PEP":
                if header is None or len(record) != len(header):
                    raise ValueError(f"Invalid mzTab peptide row at line {line_number}")
                row = dict(zip(header, record))
                if row.get("opt_global_cv_MS:1002217_decoy_peptide") == "1":
                    continue
                identity = (row["sequence"], row["modifications"])
                peptides.setdefault(identity, {
                    "sequence": identity[0], "modifications": identity[1],
                    "peptidoform": row.get("opt_global_cv_MS:1000889_peptidoform_sequence", "null"),
                })
                for key, observation in read_ratio_arrays(row, PEPTIDE_PREFIX, "_channel", normalized):
                    add_ratio(ratios, (identity, *key), observation)
    if header is None:
        raise ValueError("No mzTab peptide section found")
    return peptides, ratios


def read_proteins(root, normalized):
    proteins, ratios = {}, {}
    for identification in root.findall(".//ProteinIdentification"):
        params = user_params(identification)
        hits = {hit.get("id"): hit.get("accession") for hit in identification.findall("ProteinHit")}
        for name, value in params.items():
            if not re.fullmatch(r"indistinguishable_proteins_\d+", name):
                continue
            members = value.split(",")[1:]
            if not members or any(member not in hits for member in members):
                raise ValueError(f"Unresolved protein group members: {name}")
            identity = tuple(sorted({hits[member] for member in members}))
            proteins[identity] = {"protein_group": ";".join(identity)}
            for key, observation in read_ratio_arrays(
                params, name + "_fraction_group_level_ratio", "_label", normalized
            ):
                add_ratio(ratios, (identity, *key), observation)
    if not proteins:
        raise ValueError("No inferred protein groups found in consensusXML")
    return proteins, ratios


def make_rows(entities, ratios, groups, samples, reference, labels, normalized):
    for _, group_id, channel, _ in ratios:
        if group_id not in groups or channel not in groups[group_id]["samples"]:
            raise ValueError(f"Ratio group/channel {group_id}/{channel} is absent from the design")
    for identity, entity in sorted(entities.items()):
        for group_id, group in sorted(groups.items()):
            for channel, sample in sorted(group["samples"].items()):
                if channel == reference:
                    continue
                for normalization in (["raw", "normalized"] if normalized else ["raw"]):
                    value, count = ratios.get((identity, group_id, channel, normalization), (None, None))
                    status = "reported"
                    log2_ratio = None
                    if value is None:
                        status = "not_reported"
                    elif not math.isfinite(value):
                        status = "non_finite"
                        value = None
                    elif value <= 0:
                        status = "non_positive"
                    else:
                        log2_ratio = math.log2(value)
                    reference_sample = group["samples"][reference]
                    yield entity | {
                        "fraction_group": group_id,
                        "runs": ";".join(dict.fromkeys(run for _, run in sorted(group["runs"]))),
                        "numerator_channel": channel, "reference_channel": reference,
                        "numerator_label_definition": f"[{labels[channel - 1]}]",
                        "reference_label_definition": f"[{labels[reference - 1]}]",
                        "numerator_sample": sample, "reference_sample": reference_sample,
                        "numerator_condition": samples[sample].get("MSstats_Condition", "null"),
                        "reference_condition": samples[reference_sample].get("MSstats_Condition", "null"),
                        "normalization": normalization, "ratio": value, "log2_ratio": log2_ratio,
                        "ratio_count": count, "ratio_status": status,
                    }


def export_ratios(mztab, consensusxml, design):
    root = ET.parse(consensusxml).getroot()
    reference, labels, normalized = read_settings(root)
    groups, samples = read_design(design)
    for group_id, group in groups.items():
        if reference not in group["samples"] or any(
            not 1 <= channel <= len(labels) for channel in group["samples"]
        ):
            raise ValueError(f"Design channels do not match executed labels/reference: group {group_id}")
    check_design_maps(root, groups)
    peptide_entities, peptide_ratios = read_peptides(mztab, normalized)
    protein_entities, protein_ratios = read_proteins(root, normalized)
    # Resolve and validate both exports before writing either table.
    arguments = (groups, samples, reference, labels, normalized)
    return (list(make_rows(peptide_entities, peptide_ratios, *arguments)),
            list(make_rows(protein_entities, protein_ratios, *arguments)))


def write_table(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: "null" if value is None else value for key, value in row.items()}
                         for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mztab", type=Path, required=True)
    parser.add_argument("--consensusxml", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        peptides, proteins = export_ratios(args.mztab, args.consensusxml, args.design)
    except (ValueError, KeyError, OSError, ET.ParseError) as error:
        parser.exit(1, f"Cannot export SILAC ratios: {error}\n")
    args.outdir.mkdir(parents=True, exist_ok=True)
    for filename, fields, rows in (
        ("peptide_ratios.tsv", PEPTIDE_FIELDS, peptides),
        ("protein_ratios.tsv", PROTEIN_FIELDS, proteins),
    ):
        write_table(args.outdir / filename, fields, rows)
        reported = sum(row["ratio_status"] == "reported" for row in rows)
        print(f"{args.outdir / filename}: {reported}/{len(rows)} ratios reported")


if __name__ == "__main__":
    main()
