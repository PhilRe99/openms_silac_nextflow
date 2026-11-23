process SAVE_SILAC_CHEMISTRY {

    input:
    path chemistry

    output:
    path "metadata/silac_chemistry.json", emit: chemistry

    script:
    """
    mkdir -p metadata
    cp "${chemistry}" metadata/silac_chemistry.json
    """
    
}