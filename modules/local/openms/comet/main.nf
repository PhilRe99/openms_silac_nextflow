process COMET {

    tag "${meta.id}"

    container params.openms_container

    input: 
    tuple val(meta), path(mzml), path(database)

    output:
    tuple val(meta), path(mzml), path("${meta.id}.comet.idXML"), emit: runs

    script:

    def args = task.ext.args ?: ''

    def fragmentTolerance = params.comet_fragment_mass_tolerance != null
        ? params.comet_fragment_mass_tolerance
        : meta.fragment_mass_tolerance

    def fragmentUnit = params.comet_fragment_mass_tolerance != null
        ? 'Da'
        : meta.fragment_mass_tolerance_unit
        
    if (fragmentUnit != 'Da') {
        error "Comet requires fragment tolerance in Da; ${meta.id} specifies ${fragmentTolerance} ${fragmentUnit}. Set --comet_fragment_mass_tolerance explicitly in Da (ppm has no fixed Da equivalent)."
    }
    if (!(fragmentTolerance.toString().isNumber()) || fragmentTolerance.toString().toDouble() < 0.0001) {
        error "Comet fragment tolerance must be a number >= 0.0001 Da; got ${fragmentTolerance}."
    }

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
        ${args} \
        -out "${meta.id}.comet.idXML" \
        -comet_executable comet.exe \
        -enzyme '${meta.enzyme.toString().replace("'", "'\"'\"'")}' \
        -fragment_mass_tolerance '${fragmentTolerance}' \
        -fragment_error_units 'Da' \
        -precursor_mass_tolerance '${meta.precursor_mass_tolerance.toString().replace("'", "'\"'\"'")}' \
        -precursor_error_units '${meta.precursor_mass_tolerance_unit.toString().replace("'", "'\"'\"'")}' \
        -fixed_modifications '${meta.fixed_modifications.toString().replace("'", "'\"'\"'")}' \
        -variable_modifications "\${VARIABLE_MODS[@]}" \
        -binary_modifications "\${BINARY_MODS[@]}" \
        -reindex false \
        -threads ${task.cpus} \
        -force
    """

}
