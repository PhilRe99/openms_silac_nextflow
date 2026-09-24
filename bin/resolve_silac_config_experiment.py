import csv
import json
import sys
from pathlib import Path

CHANNEL_CHEMISTRY = {
    "light": [],
    "medium": ["Lys4", "Arg6"],
    "heavy": ["Lys8", "Arg10"]
}

#variable modification format required for comet
LABEL_MODIFICATIONS = {
    "Lys4": "Label:2H(4) (K)",
    "Arg6": "Label:13C(6) (R)",
    "Lys8": "Label:13C(6)15N(2) (K)",
    "Arg10": "Label:13C(6)15N(4) (R)",
}


#reads openms.tsv and selects current run to create json
def read_openms_data(openms_data_file, raw_data_file) -> dict:

    filename = Path(raw_data_file).name
    design = None

    with open(openms_data_file, "r") as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            if row["URI"] == filename and row["Label"] == "SILAC":
                design = row

    if design is None:
        raise ValueError(f"No SILAC OpenMS row found for {filename}")

    return design

#reads experimental_design.tsv to find row corresponding to current run rows and create json
def read_experimental_design(experimental_design_file, raw_data_file):

    with open(experimental_design_file, "r") as f:
        lines = list(f)

    matching_rows = []
    filename = Path(raw_data_file).name

    #finding seperator for the two tables.
    separator = next(i for i, line in enumerate(lines)if not line.strip())

    file_reader = csv.DictReader(lines[:separator], delimiter="\t")
    sample_reader = csv.DictReader(lines[separator + 1:],delimiter="\t")

    file_rows = list(file_reader)
    sample_rows = list(sample_reader)

    for row in file_rows:
        if Path(row["Spectra_Filepath"]).name == filename:
                matching_rows.append(row)

    if not matching_rows:
        raise ValueError(f"No matching row found for {filename} in experimental design file {experimental_design_file}")

    return matching_rows, sample_rows


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

#determines modification from label. e.g: 1,3 = [], ["light", "heavy"]
#then gets channel_chemistry e.g.: medium = ["Lys4", "Arg6"]
#then gets comet readable modifications Lys8 = Label:13C(6)15N(2) (K)
#generates json
def map_labels_to_modifications(labels, label_channels) -> dict:
    label_modifications = {}

    for label in labels:

        channel_name = label_channels[label]
        label_modifications[channel_name] = []

        for modification in CHANNEL_CHEMISTRY[channel_name]:
            label_modifications[channel_name].append(LABEL_MODIFICATIONS[modification])

    return label_modifications

#builds Comet variable and binary modifications
#generic variable mods use group 0
#each labelled SILAC channel gets its own non-zero binary group
def build_comet_binary_modifications(openms_design, label_modifications) -> list:

    variable_modifications = openms_design.get("VariableModifications", [])
    variable_modifications_silac = []
    binary_modifications = []

    for key in label_modifications.keys():
        variable_modifications_silac.extend(label_modifications[key])

    binary_group = 0

    for entry in variable_modifications.split(","):
        if entry == "": break
        binary_modifications.append(binary_group)

    for entry in label_modifications.values():
        if entry == []:
                continue
        binary_group += 1
        for value in entry:
            binary_modifications.append(binary_group)
        
    variable_modifications = variable_modifications.split(",") + variable_modifications_silac

    return variable_modifications, binary_modifications

#builds labels for FeatureFinderMultiplex in needed format
def build_ffm_labels(labels, label_channels) -> str:

    channel_names = [label_channels[label]for label in labels]

    channel_chemistries = [CHANNEL_CHEMISTRY[name]for name in channel_names]

    ffm_labels = ""

    for chem in channel_chemistries:
        ffm_labels += "[" + ",".join(chem) + "]"

    return ffm_labels

def main():

    openms_data = sys.argv[1]
    experimental_design = sys.argv[2]
    raw_data = sys.argv[3]
    silac_config = sys.argv[4]


    openms_design = read_openms_data(openms_data, raw_data)
    experimental_design_rows, _ = read_experimental_design(experimental_design, raw_data)

    labels = collect_labels(experimental_design_rows)

    global_labels = collect_global_labels(experimental_design)
    global_label_channels = determine_label_channels(global_labels)

    label_channels = {label: global_label_channels[label] for label in labels}

    label_modifications = map_labels_to_modifications(labels, label_channels)
    variable_modifications, binary_modifications = build_comet_binary_modifications(openms_design, label_modifications)
    ffm_labels = build_ffm_labels(labels, label_channels)
    
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
        "variable_modifications": variable_modifications,
        "binary_modifications": binary_modifications,
        "ffm_labels": ffm_labels
    }

    with open(silac_config, "w") as f:
        json.dump(config, f, indent=4)

if __name__ == "__main__":
    main()