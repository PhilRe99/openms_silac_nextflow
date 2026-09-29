process COMET {

    tag "${meta.id}"

    input: 
    tuple val(meta), path(mzml), path(config), path(database)

    output:
    tuple val(meta), path(mzml), path("${meta.id}.comet.idXML"), path(config), emit: runs

    script:
    """
    FIXED_MODIFICATIONS="\$(jq -r '.openms_design.FixedModifications' "${config}")"
    ENZYME="\$(jq -r '.openms_design.Enzyme' "${config}")"

    PRECURSOR_TOLERANCE="\$(jq -r '.openms_design.PrecursorMassTolerance' "${config}")"
    PRECURSOR_UNIT="\$(jq -r '.openms_design.PrecursorMassToleranceUnit' "${config}")"

    FRAGMENT_TOLERANCE="\$(jq -r '.openms_design.FragmentMassTolerance' "${config}")"
    FRAGMENT_UNIT="\$(jq -r '.openms_design.FragmentMassToleranceUnit' "${config}")"

    mapfile -t VARIABLE_MODS < <(
        jq -r '.variable_modifications[]' "${config}"
    )

    mapfile -t BINARY_MODS < <(
        jq -r '.binary_modifications[]' "${config}"
    )

    CometAdapter \
        -in "${mzml}" \
        -database "${database}" \
        -out "${meta.id}.comet.idXML" \
        -comet_executable /usr/share/OpenMS/THIRDPARTY/Comet/comet.exe \
        -enzyme "\${ENZYME}" \
        -missed_cleavages 2 \
        -fragment_mass_tolerance "\${FRAGMENT_TOLERANCE}" \
        -fragment_error_units "\${FRAGMENT_UNIT}" \
        -precursor_mass_tolerance "\${PRECURSOR_TOLERANCE}" \
        -precursor_error_units "\${PRECURSOR_UNIT}" \
        -instrument high_res \
        -fixed_modifications "\${FIXED_MODIFICATIONS}" \
        -variable_modifications "\${VARIABLE_MODS[@]}" \
        -binary_modifications "\${BINARY_MODS[@]}" \
        -isotope_error 0/1/2/3 \
        -reindex false \
        -threads ${task.cpus} \
        -force
    """

}

