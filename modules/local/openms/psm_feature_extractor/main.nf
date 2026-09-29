process PSM_FEATURE_EXTRACTOR {

    tag "${meta.id}"

    input:
    tuple val(meta), path(mzml), path(idxml), path(config)

    output:
    tuple val(meta), path(mzml), 
        path("${meta.id}.comet.features.idXML"),
        path(config), emit: runs

    
    script:
    """
    PSMFeatureExtractor \
        -in "${idxml}" \
        -out "${meta.id}.comet.features.idXML"
    """

}
