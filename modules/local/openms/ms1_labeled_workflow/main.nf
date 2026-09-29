process MS1_LABELED_WORKFLOW {

    input:
    tuple val(labels), path(mzmls), path(idxmls)
    path design
    path database

    output:
    path "experiment.mzTab", emit: mztab
    path "experiment.consensusXML", emit: consensus
    

    script:

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
    -out experiment.mzTab \
    -out_cxml experiment.consensusXML \
    -threads ${task.cpus}
    """
}