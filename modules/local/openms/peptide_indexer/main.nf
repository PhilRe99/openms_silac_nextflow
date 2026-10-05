process PEPTIDE_INDEXER {

    tag "${meta.id}"

    container params.openms_container
    
    input:
    tuple val(meta), path(mzml), path(idxml)
    path(database)

    output:
    tuple val(meta), path(mzml), path("${meta.id}.comet.indexed.idXML"), emit: runs
    
    script:
    """
    PeptideIndexer \
        -in "${idxml}" \
        -fasta "${database}" \
        -out "${meta.id}.comet.indexed.idXML" \
        -decoy_string DECOY_ \
        -decoy_string_position prefix \
        -threads ${task.cpus}
    """
}
