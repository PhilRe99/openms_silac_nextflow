process GENERATE_DECOY_DATABASE {

    tag "${fasta.baseName}"

    input:
    path fasta

    output:
    path "${fasta.baseName}_target_decoy.fasta", emit: database

    script:
    """
    DecoyDatabase \
        -in "${fasta}" \
        -out "${fasta.baseName}_target_decoy.fasta" \
        -decoy_string DECOY_ \
        -decoy_string_position prefix \
        -method reverse
    """

}
