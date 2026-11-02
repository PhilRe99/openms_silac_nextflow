process SAVE_SILAC_CHEMISTRY {

    container 'openms-silac/python:3.11-procps'

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