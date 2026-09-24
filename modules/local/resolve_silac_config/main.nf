process RESOLVE_SILAC_CONFIG {

    tag "${meta.id}"

    input:
    tuple val(meta), path(mzml)
    path openms_config
    path experimental_design

    output:
    tuple val(meta), path(mzml), path("silac_config_${meta.id}.json"), emit: runs

    script:
    """
    python3 "${projectDir}/bin/resolve_silac_config_experiment.py" \
        "${openms_config}" \
        "${experimental_design}" \
        "${mzml}" \
        "silac_config_${meta.id}.json"
    """


}