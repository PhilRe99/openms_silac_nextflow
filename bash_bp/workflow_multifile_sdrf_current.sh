#!/usr/bin/env bash
set -euo pipefail

#Input args
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <mzML_dir> <fasta> <sdrf> <execution_id>" >&2
    exit 1
fi

INPUT_DIR="$1"
FASTA="$2"
SDRF_METADATA="$3"
EXECUTION_ID="$4"

#Validate args present and filetypes
[[ -d "$INPUT_DIR" ]] || { echo "ERROR: mzML directory not found: $INPUT_DIR" >&2; exit 1; }
[[ -f "$FASTA" ]] || { echo "ERROR: FASTA not found: $FASTA" >&2; exit 1; }
[[ -f "$SDRF_METADATA" ]] || { echo "ERROR: SDRF not found: $SDRF_METADATA" >&2; exit 1; }


#Derive names for workflow outputs for dynamic naming
SDRF_METADATA="$(realpath "$SDRF_METADATA")"

DB_FILENAME="$(basename "$FASTA")"
DB_BASENAME="${DB_FILENAME%.*}"
METADATA_FILENAME="$(basename "$SDRF_METADATA")"
EXPERIMENT_ID="${METADATA_FILENAME%.sdrf.tsv}"

RUN_NAMESPACE="${EXPERIMENT_ID}/current/${EXECUTION_ID}"
WORK_DIR="work/${RUN_NAMESPACE}"
RESULTS_DIR="results/${RUN_NAMESPACE}"
LOG_DIR="logs/${RUN_NAMESPACE}"
METADATA_DIR="${WORK_DIR}/metadata"
CENTROIDED_DIR="${WORK_DIR}/centroided"


mkdir -p "$WORK_DIR" "$METADATA_DIR" "$RESULTS_DIR" "$CENTROIDED_DIR" "$LOG_DIR"


#current hardcoded vals
MISSED_CLEAVAGES=2
PSM_FDR=0.01
PROTEIN_FDR=0.01

#out experimental design.tsv and openms.tsv
parse_sdrf validate-sdrf   --sdrf_file "$SDRF_METADATA"   --template ms-proteomics

(
    cd "$METADATA_DIR"
    parse_sdrf convert-openms -s "$SDRF_METADATA"
)

#parse_sdrf convert-openms   -s $SDRF_METADATA

EXPERIMENTAL_DESIGN="${METADATA_DIR}/experimental_design.tsv"
OPENMS_TSV="${METADATA_DIR}/openms.tsv"
TARGET_DECOY_DB="${WORK_DIR}/${DB_BASENAME}_target_decoy.fasta"

DecoyDatabase \
  -in "$FASTA" \
  -out "$TARGET_DECOY_DB" \
  -decoy_string DECOY_ \
  -decoy_string_position prefix \
  -method reverse

readarray -t DESIGN_MZMLS < <(
    python3 python_scripts/build_input_file_array.py "$EXPERIMENTAL_DESIGN"
)


FILE_STATE="${METADATA_DIR}/files.tsv"
printf 'base\tmzml\tidxml\tconfig\tlabels\n' > "$FILE_STATE"

ALL_MZMLS=()
ALL_IDS=()
MS1_LABELS="[][Lys4,Arg6][Lys8,Arg10]"

for DESIGN_MZML in "${DESIGN_MZMLS[@]}"; do

    RAW_FILENAME="$(basename "$DESIGN_MZML")"
    RAW_BASENAME="${RAW_FILENAME%.mzML}"
    INPUT_MZML="${INPUT_DIR}/${RAW_FILENAME}"
    CENTROIDED_MZML="${CENTROIDED_DIR}/${RAW_FILENAME}"
    SILAC_CONFIG="${METADATA_DIR}/silac_config_${RAW_BASENAME}.json"

    COMET_ID="${WORK_DIR}/${RAW_BASENAME}.comet.idXML"
    INDEXED_ID="${WORK_DIR}/${RAW_BASENAME}.comet.indexed.idXML"
    FEATURES_ID="${WORK_DIR}/${RAW_BASENAME}.comet.features.idXML"
    PERCOLATOR_ID="${WORK_DIR}/${RAW_BASENAME}.percolator.pep.fdr01.idXML"

    [[ -f "$INPUT_MZML" ]] || {
        echo "ERROR: mzML from experimental design not found: $INPUT_MZML" >&2
        exit 1
    }

    echo
    echo "========== $RAW_FILENAME =========="

    #Detect profile vs centroided data w/ FileInfo
    #If profile, run PeakPickerHiRes to centroid the data if not skip

    FILEINFO_OUTPUT="$(FileInfo -in "$INPUT_MZML")"
    if grep -q "Profile (Profile)" <<< "$FILEINFO_OUTPUT"; then
        PeakPickerHiRes -in "$INPUT_MZML" -out "$CENTROIDED_MZML"
        ANALYSIS_MZML="$CENTROIDED_MZML"
    elif grep -q "Centroid (Centroid)" <<< "$FILEINFO_OUTPUT"; then
        ANALYSIS_MZML="$INPUT_MZML"
    else
        echo "ERROR: Could not determine spectrum type from FileInfo: $INPUT_MZML" >&2
        exit 1
    fi

    python3 python_scripts/resolve_silac_config_experiment.py \
        "$OPENMS_TSV" \
        "$EXPERIMENTAL_DESIGN" \
        "$INPUT_MZML" \
        "$SILAC_CONFIG"

    #assigning values directly from json file using jq for easy access
    FIXED_MODIFICATIONS="$(jq -r '.openms_design.FixedModifications' "$SILAC_CONFIG")"
    ENZYME="$(jq -r '.openms_design.Enzyme' "$SILAC_CONFIG")"

    #some of these will be defaults but still added so if they should change in any future run
    #they are correctly assigned
    PRECURSOR_TOLERANCE="$(jq -r '.openms_design.PrecursorMassTolerance' "$SILAC_CONFIG")"
    PRECURSOR_TOLERANCE_UNIT="$(jq -r '.openms_design.PrecursorMassToleranceUnit' "$SILAC_CONFIG")"

    FRAGMENT_TOLERANCE="$(jq -r '.openms_design.FragmentMassTolerance' "$SILAC_CONFIG")"
    FRAGMENT_TOLERANCE_UNIT="$(jq -r '.openms_design.FragmentMassToleranceUnit' "$SILAC_CONFIG")"

    FFM_LABELS="$(jq -r '.ffm_labels' "$SILAC_CONFIG")"

    #mapfile: reads lines of input and stores each line as one element of a Bash array
    #variable_modifications[] prints one element per line
    #-t removes \n
    mapfile -t COMET_VARIABLE_MODS < <(
        jq -r '.variable_modifications[]' "$SILAC_CONFIG"
    )

    mapfile -t COMET_BINARY_MODS < <(
        jq -r '.binary_modifications[]' "$SILAC_CONFIG"
    )

    CometAdapter \
        -in "$ANALYSIS_MZML" \
        -database "$TARGET_DECOY_DB" \
        -out "$COMET_ID" \
        -comet_executable /usr/share/OpenMS/THIRDPARTY/Comet/comet.exe \
        -enzyme "$ENZYME" \
        -missed_cleavages "$MISSED_CLEAVAGES" \
        -fragment_mass_tolerance "$FRAGMENT_TOLERANCE" \
        -fragment_error_units "$FRAGMENT_TOLERANCE_UNIT" \
        -precursor_mass_tolerance "$PRECURSOR_TOLERANCE" \
        -precursor_error_units "$PRECURSOR_TOLERANCE_UNIT" \
        -instrument high_res \
        -fixed_modifications "$FIXED_MODIFICATIONS" \
        -variable_modifications "${COMET_VARIABLE_MODS[@]}" \
        -binary_modifications "${COMET_BINARY_MODS[@]}" \
        -isotope_error 0/1/2/3 \
        -reindex false \
        -threads 0 \
        -force

    PeptideIndexer \
        -in "$COMET_ID" \
        -fasta "$TARGET_DECOY_DB" \
        -out "$INDEXED_ID" \
        -decoy_string DECOY_ \
        -decoy_string_position prefix

    #optional can also just use flag at indexer
    PSMFeatureExtractor \
        -in "$INDEXED_ID" \
        -out "$FEATURES_ID"

    PercolatorAdapter \
        -in "$FEATURES_ID" \
        -out "$PERCOLATOR_ID" \
        -use_subprocess true \
        -percolator_executable /usr/share/OpenMS/THIRDPARTY/Percolator/percolator \
        -score_type pep \
        -score:fdr "$PSM_FDR"

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$RAW_BASENAME" "$ANALYSIS_MZML" "$PERCOLATOR_ID" "$SILAC_CONFIG" "$FFM_LABELS" \
        >> "$FILE_STATE"

    ALL_MZMLS+=("$ANALYSIS_MZML")
    ALL_IDS+=("$PERCOLATOR_ID")

done

echo "mzML inputs:"
printf '  %s\n' "${ALL_MZMLS[@]}"

echo "ID inputs:"
printf '  %s\n' "${ALL_IDS[@]}"

echo "Labels: $MS1_LABELS"
echo

MS1LabeledWorkflow \
    -in "${ALL_MZMLS[@]}" \
    -ids "${ALL_IDS[@]}" \
    -labels "[][Lys4,Arg6][Lys8,Arg10]" \
    -design "$EXPERIMENTAL_DESIGN" \
    -fasta "$TARGET_DECOY_DB" \
    -out "$RESULTS_DIR/experiment.mzTab" \
    -out_cxml "$RESULTS_DIR/experiment.consensusXML" \
    -threads 0

#    -max_nr_labelled_aas "$MISSED_CLEAVAGES" \
#    -proteinFDR "$PROTEIN_FDR" \
#    -picked_proteinFDR true \
#    -ProteinQuantification:fractions:aggregate sum \
#    -ratios:reference_channel 1 \
#    -ratios:min_ratio_count 2 \
#    -ratios:normalize false \

echo "done"
