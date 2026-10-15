import argparse
import csv
import json
import math
import re
from pathlib import Path

LABEL_MODIFICATIONS = {
    "Lys4": "Label:2H(4) (K)",
    "Arg6": "Label:13C(6) (R)",
    "Lys6": "Label:13C(6) (K)",   
    "Lys8": "Label:13C(6)15N(2) (K)",
    "Arg10": "Label:13C(6)15N(4) (R)",
}
SILAC_LABEL_MODIFICATIONS = set(LABEL_MODIFICATIONS.values())

#channel numbers written by sdrf-pipelines into the experimental design
#(OpenMSExperimentalDesign._get_label_from_labels): light/medium/heavy = 1/2/3 as soon as
#any SDRF row is 'SILAC medium', otherwise light/heavy = 1/2
SDRF_SILAC_3PLEX = {"silac light": "1", "silac medium": "2", "silac heavy": "3"}
SDRF_SILAC_2PLEX = {"silac light": "1", "silac heavy": "2"}

def split_modifications(value) -> list[str]:
    return [mod.strip() for mod in (value or "").split(",") if mod.strip()]

def run_stem(filename):
    return re.sub(r"(?i)\.(raw|mzml)$", "", Path(filename).name)

def resolve_comet_fragment_tolerance(metadata_value, metadata_unit, override=None):
    """Resolve the Comet Da bin tolerance without changing the input metadata."""
    value = override if override is not None else metadata_value
    unit = "Da" if override is not None else metadata_unit
    source = "explicit_override" if override is not None else "metadata"

    if unit != "Da":
        raise ValueError(
            f"Comet requires fragment tolerance in Da; metadata specifies "
            f"{metadata_value} {metadata_unit}. Set --comet_fragment_mass_tolerance "
            "explicitly in Da (ppm has no fixed Da equivalent)."
        )
    try:
        tolerance = float(value)
    except (TypeError, ValueError):
        raise ValueError(
            f"Comet fragment tolerance must be a finite number >= 0.0001 Da; got {value!r}."
        ) from None
    if isinstance(value, bool) or not math.isfinite(tolerance) or tolerance < 0.0001:
        raise ValueError(
            f"Comet fragment tolerance must be a finite number >= 0.0001 Da; got {value!r}."
        )

    return {
        "value": tolerance,
        "unit": unit,
        "source": source,
        "metadata_value": metadata_value,
        "metadata_unit": metadata_unit,
    }

#reads openms.tsv and selects current run to create json
def read_openms_data(openms_data_file, raw_data_file) -> dict:

    filename = Path(raw_data_file).name
    design = None

    with open(openms_data_file, "r") as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            if row["Filename"] == filename and row["Label"] == "SILAC":
                design = row

    if design is None:
        raise ValueError(f"No SILAC OpenMS row found for {filename}")

    return design

#reads experimental_design.tsv to find row corresponding to current run rows and create json
def read_experimental_design(experimental_design_file, raw_data_file):

    with open(experimental_design_file, "r") as f:
        lines = list(f)

    filename = Path(raw_data_file).name
    design_filenames = [filename]
    if Path(filename).suffix.lower() == ".raw":
        # The converter can write mzML design paths before RAW conversion runs.
        design_filenames.append(Path(filename).with_suffix(".mzML").name)

    #finding seperator for the two tables.
    separator = next(i for i, line in enumerate(lines)if not line.strip())

    file_reader = csv.DictReader(lines[:separator], delimiter="\t")
    file_rows = list(file_reader)


    for design_filename in design_filenames:
        matching_rows = [
            row for row in file_rows
            if Path(row["Spectra_Filepath"]).name == design_filename
        ]
        if matching_rows:
            return matching_rows

    raise ValueError(f"No matching row found for {filename} in experimental design file {experimental_design_file}")


#finds needed labels. e.g.: light, heavy = 1,3
def collect_labels(experimental_design_rows) -> list[str]:
    
    labels = set()

    for row in experimental_design_rows:
        labels.add(row["Label"])

    return sorted(labels, key=int)

#finds global labels
#triplex 1,2,3
#duplex 1,2 
#so if any triplex is included it generates 1,2,3 which then leads to all duplex files using 1,3
#if only duplex runs it is 1,2
def collect_global_labels(experimental_design_file):

    with open(experimental_design_file, "r") as f:
        lines = list(f)

    separator = next(
        i for i, line in enumerate(lines)
        if not line.strip()
    )

    file_reader = csv.DictReader(
        lines[:separator],
        delimiter="\t"
    )

    labels = {
        row["Label"]
        for row in file_reader
    }

    return labels

#checks global labels to determine what label is finally chosen
def determine_label_channels(global_labels):

    if global_labels == {"1", "2"}:
        return {
            "1": "light",
            "2": "heavy"
        }

    elif global_labels == {"1", "2", "3"}:
        return {
            "1": "light",
            "2": "medium",
            "3": "heavy"
        }

    else:
        raise ValueError(
            f"Unsupported SILAC label configuration: "
            f"{sorted(global_labels, key=int)}"
        )

#reads the channel names of the current run from the SDRF, e.g.: {"1": "light", "2": "medium"}
#the design numbers alone cannot tell a light/medium from a light/heavy duplex
def read_sdrf_label_channels(sdrf_file, raw_data_file) -> dict:

    with open(sdrf_file, newline="") as f:
        rows = [row for row in csv.reader(f, delimiter="\t") if row]

    header = [column.strip().lower() for column in rows[0]]
    for column in ("comment[data file]", "comment[label]"):
        if column not in header:
            raise ValueError(f"SDRF {sdrf_file} has no {column} column")
    file_column = header.index("comment[data file]")
    label_column = header.index("comment[label]")

    records = []
    for row in rows[1:]:
        #plain ('SILAC heavy') or NT-encoded ('NT=SILAC heavy;AC=...') label
        match = re.search(r"(?:^|;)\s*NT=([^;]+)", row[label_column], flags=re.IGNORECASE)
        label = (match.group(1) if match else row[label_column]).strip().lower()
        records.append((row[file_column].strip(), label))

    silac_labels = {label for _, label in records if label.startswith("silac")}
    numbering = SDRF_SILAC_3PLEX if "silac medium" in silac_labels else SDRF_SILAC_2PLEX

    filename = Path(raw_data_file).name
    run_labels = [label for data_file, label in records if data_file == filename]
    if not run_labels:
        run_labels = [label for data_file, label in records if run_stem(data_file) == run_stem(filename)]
    if not run_labels:
        raise ValueError(f"No SDRF row found for {filename} in {sdrf_file}")

    channels = {}
    for label in run_labels:
        if label not in numbering:
            raise ValueError(f"Unsupported SILAC label {label!r} for {filename} in {sdrf_file}")
        channels[numbering[label]] = label.removeprefix("silac ")

    return channels

#determines modification from label. e.g: 1,3 = [], ["light", "heavy"]
#then gets channel_chemistry e.g.: medium = ["Lys4", "Arg6"]
#then gets comet readable modifications Lys8 = Label:13C(6)15N(2) (K)
#generates json
def map_labels_to_modifications(labels, label_channels, channel_chemistry) -> dict:
    label_modifications = {}

    for label in labels:

        channel_name = label_channels[label]
        if channel_name not in channel_chemistry:
            raise ValueError(f"SILAC channel {channel_name!r} is not defined in the chemistry file")
        label_modifications[channel_name] = []

        for modification in channel_chemistry[channel_name]:

            if modification not in LABEL_MODIFICATIONS:
                raise ValueError(
                    f"Unsupported SILAC label {modification!r} "
                    f"in channel {channel_name!r}. "
                    f"Supported labels: {', '.join(LABEL_MODIFICATIONS)}"
                )

            label_modifications[channel_name].append(LABEL_MODIFICATIONS[modification])

    return label_modifications

#builds Comet variable and binary modifications
#generic variable mods use group 0
#each labelled SILAC channel gets its own non-zero binary group
#SILAC labels the SDRF lists as variable modifications are searched in their binary group only:
#an additional group 0 copy would also allow partially labelled peptides
def build_comet_binary_modifications(openms_design, label_modifications) -> list:

    sdrf_modifications = split_modifications(openms_design.get("VariableModifications"))

    variable_modifications_silac = []

    for key in label_modifications.keys():
        variable_modifications_silac.extend(label_modifications[key])

    declared_labels = {mod for mod in sdrf_modifications if mod in SILAC_LABEL_MODIFICATIONS}
    if declared_labels and declared_labels != set(variable_modifications_silac):
        raise ValueError(
            f"The SDRF declares the SILAC labels {sorted(declared_labels)}, but the chemistry file "
            f"resolves {sorted(set(variable_modifications_silac))} for this run. "
            "Correct the SDRF modification parameters or --silac_chemistry."
        )

    variable_modifications = [mod for mod in sdrf_modifications if mod not in SILAC_LABEL_MODIFICATIONS]

    binary_group = 0

    binary_modifications = [0] * len(variable_modifications)

    for entry in label_modifications.values():
        if entry == []:
                continue
        binary_group += 1
        for value in entry:
            binary_modifications.append(binary_group)
        
    variable_modifications = variable_modifications + variable_modifications_silac
    
    return variable_modifications, binary_modifications

#builds Comet fixed modifications as a list, one entry per modification
def build_fixed_modifications(openms_design) -> list:

    fixed_modifications = split_modifications(openms_design.get("FixedModifications"))

    fixed_labels = [mod for mod in fixed_modifications if mod in SILAC_LABEL_MODIFICATIONS]
    if fixed_labels:
        raise ValueError(
            f"The SDRF annotates the SILAC labels {fixed_labels} as fixed modifications, "
            "which would label every channel. Annotate them as variable (MT=variable)."
        )

    return fixed_modifications

#builds labels for FeatureFinderMultiplex in needed format
def build_ffm_labels(labels, label_channels, channel_chemistry) -> str:

    channel_names = [label_channels[label] for label in labels]
    channel_chemistries = [channel_chemistry[name] for name in channel_names]

    ffm_labels = ""

    for chem in channel_chemistries:
        ffm_labels += "[" + ",".join(chem) + "]"

    return ffm_labels

def main(argv=None):
    parser = argparse.ArgumentParser(description="Resolve per-run SILAC and Comet fragment settings.")
    parser.add_argument("openms_data")
    parser.add_argument("experimental_design")
    parser.add_argument("raw_data")
    parser.add_argument("silac_config")
    parser.add_argument("chemistry_file")
    parser.add_argument("--comet-fragment-mass-tolerance", default=None,
                        help="Explicit Comet fragment bin tolerance override, always in Da.")
    parser.add_argument("--sdrf", default=None,
                        help="SDRF the design was converted from. Its comment[label] names the channels; "
                             "without it a light/medium duplex is taken for light/heavy.")
    args = parser.parse_args(argv)

    openms_data = args.openms_data
    experimental_design = args.experimental_design
    raw_data = args.raw_data
    silac_config = args.silac_config
    chemistry_file = args.chemistry_file

    
    with open(chemistry_file) as handle:
        channel_chemistry = json.load(handle)

    openms_design = read_openms_data(openms_data, raw_data)
    
    try:
        fragment_tolerance = resolve_comet_fragment_tolerance(
            openms_design.get("FragmentMassTolerance"),
            openms_design.get("FragmentMassToleranceUnit"),
            args.comet_fragment_mass_tolerance,
        )
    except ValueError as exc:
        parser.error(f"{Path(raw_data).name}: {exc}")
    experimental_design_rows = read_experimental_design(experimental_design, raw_data)

    labels = collect_labels(experimental_design_rows)

    try:
        if args.sdrf:
            label_channels = read_sdrf_label_channels(args.sdrf, raw_data)
            if set(label_channels) != set(labels):
                raise ValueError(
                    f"SDRF channels {label_channels} do not match the design labels {labels}; "
                    "the channel numbering of sdrf-pipelines may have changed."
                )
        else:
            global_labels = collect_global_labels(experimental_design)
            global_label_channels = determine_label_channels(global_labels)

            label_channels = {label: global_label_channels[label] for label in labels}

        label_modifications = map_labels_to_modifications(labels, label_channels, channel_chemistry)

        variable_modifications, binary_modifications = build_comet_binary_modifications(openms_design, label_modifications)
        fixed_modifications = build_fixed_modifications(openms_design)
        ffm_labels = build_ffm_labels(labels, label_channels, channel_chemistry)
    except ValueError as exc:
        parser.error(f"{Path(raw_data).name}: {exc}")

    #print(f"OpenMS design: {openms_design}")
    #print(f"Experimental design rows: {experimental_design_rows}")
    #print(f"Labels: {labels}")
    #print(f"Label modifications: {label_modifications}")
    #print(f"Variable modifications: {variable_modifications}")
    #print(f"Binary modifications: {binary_modifications}")
    #print(f"FeatureFinderMultiplex labels: {ffm_labels}")

    config = {
        "openms_design": openms_design,
        "experimental_design_rows": experimental_design_rows,
        "labels": labels,
        "label_channels": label_channels,
        "label_modifications": label_modifications,
        "fixed_modifications": fixed_modifications,
        "variable_modifications": variable_modifications,
        "binary_modifications": binary_modifications,
        "ffm_labels": ffm_labels,
        "comet_fragment_mass_tolerance": fragment_tolerance,
    }

    with open(silac_config, "w") as f:
        json.dump(config, f, indent=4)

if __name__ == "__main__":
    main()
