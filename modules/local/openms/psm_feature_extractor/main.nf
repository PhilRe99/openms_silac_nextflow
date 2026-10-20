process PSM_FEATURE_EXTRACTOR {

    tag "${meta.id}"

    container params.openms_container

    input:
    tuple val(meta), path(mzml), path(idparquet)

    output:
    tuple val(meta), path(mzml), 
        path("${meta.id}.comet.features.idparquet"), emit: runs

    
    script:
    """
    PSMFeatureExtractor \
        -in "${idparquet}" \
        -out "${meta.id}.comet.features.idparquet"
    """

}
