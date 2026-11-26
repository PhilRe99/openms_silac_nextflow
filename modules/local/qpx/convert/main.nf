process QPX_CONVERT {

    input:
    path consensus
    path sdrf
    path fasta

    output:
    path 'qpx', emit: qpx

    script:

    def quote = { value ->
        "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }

    def prefix = sdrf.name.replaceFirst(/(?i)(\.sdrf)?\.tsv$/, '')
    def accession = prefix.find(/(?i)PXD\d{6}/)?.toUpperCase()
    def accessionArg = accession ? "--project-accession ${quote.call(accession)}" : ''

    """
    qpxc convert openms-consensus \
        --consensusxml ${quote.call(consensus)} \
        --sdrf-file ${quote.call(sdrf)} \
        --fasta ${quote.call(fasta)} \
        --output-folder qpx \
        --output-prefix ${quote.call(prefix)} \
        ${accessionArg}
    """
}
