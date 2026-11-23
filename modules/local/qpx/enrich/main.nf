process QPX_ENRICH {

    input:
    tuple path(native_qpx, stageAs: 'native_qpx'), path(consensus)
    path sdrf

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
    qpxc convert openms \
        --qpx-dir ${quote.call(native_qpx)} \
        --sdrf-file ${quote.call(sdrf)} \
        --consensusxml ${quote.call(consensus)} \
        --output-folder qpx \
        --output-prefix ${quote.call(prefix)} \
        ${accessionArg}
    """
}
