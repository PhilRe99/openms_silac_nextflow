process COMET {

    tag "${meta.id}"

    input: 
    tuple val(meta), path(mzml), path(database)

    output:
    tuple val(meta), path(mzml), path("${meta.id}.comet.idXML"), emit: runs

    script:
    def variableMods = meta.variable_modifications
        .collect { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        .join(' ')
    def binaryMods = meta.binary_modifications
        .collect { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        .join(' ')

    """
    VARIABLE_MODS=( ${variableMods} )
    BINARY_MODS=( ${binaryMods} )

    CometAdapter \
        -in "${mzml}" \
        -database "${database}" \
        -out "${meta.id}.comet.idXML" \
        -comet_executable /usr/share/OpenMS/THIRDPARTY/Comet/comet.exe \
        -enzyme '${meta.enzyme.toString().replace("'", "'\"'\"'")}' \
        -missed_cleavages 2 \
        -fragment_mass_tolerance '${meta.fragment_mass_tolerance.toString().replace("'", "'\"'\"'")}' \
        -fragment_error_units '${meta.fragment_mass_tolerance_unit.toString().replace("'", "'\"'\"'")}' \
        -precursor_mass_tolerance '${meta.precursor_mass_tolerance.toString().replace("'", "'\"'\"'")}' \
        -precursor_error_units '${meta.precursor_mass_tolerance_unit.toString().replace("'", "'\"'\"'")}' \
        -instrument high_res \
        -fixed_modifications '${meta.fixed_modifications.toString().replace("'", "'\"'\"'")}' \
        -variable_modifications "\${VARIABLE_MODS[@]}" \
        -binary_modifications "\${BINARY_MODS[@]}" \
        -isotope_error 0/1/2/3 \
        -reindex false \
        -threads ${task.cpus} \
        -force
    """

}

