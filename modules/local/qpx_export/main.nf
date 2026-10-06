process QPX_EXPORT {

    //this is not used in favor of native openms qpx

    container params.openms_container

    publishDir {"${params.outdir}/${prefix}/qpx"}, mode: 'copy', saveAs: { name -> name.tokenize('/').last()}

    input:
    path consensusXML
    path sdrf
    val prefix

    output:
    path "qpx_output/*", emit: dataset

    script:
    """
    qpxc convert openms-consensus \
        --consensusxml "${consensusxml}" \
        --sdrf-file "${sdrf}" \
        --output-folder qpx_output \
        --output-prefix "${prefix}" \
        --structures feature,psm,pg,run,sample \
        --no-mudata

    qpxc validate \
        --dataset-path qpx_output \
        --strict \
        > qpx_output/validation.txt 2>&1

    qpxc --version > qpx_output/converter_version.txt
    """


}