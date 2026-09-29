process PSM_FEATURE_EXTRACTOR {

    tag "${meta.id}"

    input:
    tuple val(meta), path(mzml), path(idxml)

    output:
    tuple val(meta), path(mzml), 
        path("${meta.id}.comet.features.idXML"), emit: runs

    
    script:
    """
    PSMFeatureExtractor \
        -in "${idxml}" \
        -out "${meta.id}.comet.features.idXML"
    """

}
