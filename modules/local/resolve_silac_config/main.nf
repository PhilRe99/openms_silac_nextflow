process RESOLVE_SILAC_CONFIG {

    tag "${meta.id}"

    container 'openms-silac/python:3.11-procps'

    publishDir "${params.outdir}/${file(params.input).name.replaceFirst(/\.sdrf\.tsv$/, '')}/configuration", mode: 'copy', pattern: 'silac_config_*.json'

    input:
    tuple val(meta), path(mzml)
    path openms_config
    path experimental_design
    path chemistry
    val fragment_tolerance_options

    output:
    tuple val(meta), path(mzml), path("silac_config_${meta.id}.json"), emit: runs

    script:

    def override = fragment_tolerance_options.override
    def fragmentArgs = override != null
        ? "--comet-fragment-mass-tolerance='" + override.toString().replace("'", "'\"'\"'") + "'"
        : ''

    """
    python3 "${projectDir}/bin/resolve_silac_config_experiment.py" \
        "${openms_config}" \
        "${experimental_design}" \
        "${mzml}" \
        "silac_config_${meta.id}.json" \
        "${chemistry}" \
        ${fragmentArgs}
    """
}
