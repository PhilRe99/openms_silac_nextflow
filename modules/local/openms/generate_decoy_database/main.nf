process GENERATE_DECOY_DATABASE {

    tag "${fasta.baseName}"
    
    container params.openms_container

    input:
    path fasta
    path contaminants, stageAs: 'contaminants/*'

    output:
    path "${fasta.baseName}_target_decoy.fasta", emit: database

    script:
    
    def contaminantFasta = contaminants ? "\"${contaminants}\"" : ''

    """
    DecoyDatabase \
        -in "${fasta}" ${contaminantFasta} \
        -out "${fasta.baseName}_target_decoy.fasta" \
        -decoy_string DECOY_ \
        -decoy_string_position prefix \
        -method reverse
    """

}
