process PEPTIDE_INDEXER {
    tag "${meta.id}"

    input:
    tuple val(meta), path(mzml), path(idxml), path(config)
    path(database)

    output:
    tuple val(meta), path(mzml), path("${meta.id}.comet.indexed.idXML"), path(config), emit: runs
    
    script:
    """
    PeptideIndexer \
        -in "${idxml}" \
        -fasta "${database}" \
        -out "${meta.id}.comet.indexed.idXML" \
        -decoy_string DECOY_ \
        -decoy_string_position prefix
        -threads ${task.cpus}
    """
}
