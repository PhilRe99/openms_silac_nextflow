process PERCOLATOR {

    tag "${meta.id}"

    container params.openms_container

    input:
    tuple val(meta), path(mzml), path(idxml)

    output:
    tuple val(meta),
        path(mzml),
        path("${meta.id}.percolator.idXML"),
        emit: runs

    script:
    """
    PercolatorAdapter \
        -in "${idxml}" \
        -out "${meta.id}.percolator.idXML" \
        -use_subprocess true \
        -percolator_executable percolator \
        -score_type pep \
        -score:fdr ${params.psm_fdr} \
        -threads ${task.cpus}
    """
}
