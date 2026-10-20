process COMET {

    tag "${meta.id}"

    container params.openms_container

    input: 
    tuple val(meta), path(mzml), path(database)

    output:
    tuple val(meta), path(mzml), path("${meta.id}.comet.idparquet"), emit: runs

    script:

    def args = task.ext.args ?: ''

    def fragmentTolerance = meta.comet_fragment_mass_tolerance
    def fragmentUnit = meta.comet_fragment_mass_tolerance_unit
    
    if (fragmentTolerance == null || fragmentUnit != 'Da' || meta.comet_fragment_mass_tolerance_source == null) {
        error "Missing resolved Comet fragment settings for ${meta.id}; run RESOLVE_SILAC_CONFIG first."
    }

    def fixedMods = meta.fixed_modifications
        .collect { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        .join(' ')
    def variableMods = meta.variable_modifications
        .collect { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        .join(' ')
    def binaryMods = meta.binary_modifications
        .collect { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        .join(' ')

    """
    FIXED_MODS=( ${fixedMods} )
    VARIABLE_MODS=( ${variableMods} )
    BINARY_MODS=( ${binaryMods} )

    CometAdapter \
        -in "${mzml}" \
        -database "${database}" \
        ${args} \
        -out "${meta.id}.comet.idparquet" \
        -comet_executable comet.exe \
        -enzyme '${meta.enzyme.toString().replace("'", "'\"'\"'")}' \
        -fragment_mass_tolerance '${fragmentTolerance}' \
        -fragment_error_units '${fragmentUnit}' \
        -precursor_mass_tolerance '${meta.precursor_mass_tolerance.toString().replace("'", "'\"'\"'")}' \
        -precursor_error_units '${meta.precursor_mass_tolerance_unit.toString().replace("'", "'\"'\"'")}' \
        -fixed_modifications "\${FIXED_MODS[@]}" \
        -variable_modifications "\${VARIABLE_MODS[@]}" \
        -binary_modifications "\${BINARY_MODS[@]}" \
        -reindex false \
        -threads ${task.cpus} \
        -force
    """

}
