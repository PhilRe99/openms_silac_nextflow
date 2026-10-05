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
    -max_nr_labelled_aas 2 \
    -proteinFDR 0.01 \
    -picked_proteinFDR true \
    -ProteinQuantification:fractions:aggregate sum \
    -ratios:reference_channel 1 \
    -ratios:min_ratio_count 2 \
    -ratios:normalize false \
    -out experiment.mzTab \
    -out_cxml experiment.consensusXML \
    -threads ${task.cpus}
    """
}
