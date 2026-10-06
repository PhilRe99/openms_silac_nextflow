process MS1_LABELED_WORKFLOW {

    publishDir "${params.outdir}/${file(params.input).name.replaceFirst(/\.sdrf\.tsv$/, '')}", mode: 'copy'

    container params.openms_container

    input:
    tuple val(labels), path(mzmls), path(idxmls)
    path design
    path database

    output:
    path "experiment.mzTab", emit: mztab
    path "experiment.consensusXML", emit: consensus
    path "qpx", optional: true, emit: qpx
    
    script:

    def args = task.ext.args ?: ''
    def qpxArgs = params.enable_qpx ? '-out_qpx qpx' : ''

    def mzmlArgs = mzmls.collect { item ->
        "'" + item.toString().replace("'", "'\"'\"'") + "'"
    }.join(' ')
    def idxmlArgs = idxmls.collect { item ->
        "'" + item.toString().replace("'", "'\"'\"'") + "'"
    }.join(' ')

    """
    MS1LabeledWorkflow \
        -in ${mzmlArgs} \
        -ids ${idxmlArgs} \
        -labels "${labels}" \
        -design "${design}" \
        -fasta "${database}" \
        ${args} \
        -out experiment.mzTab \
        -out_cxml experiment.consensusXML \
        ${qpxArgs} \
        -threads ${task.cpus}
    """
}
