process PERCOLATOR {

    tag "${meta.id}"

    container params.openms_container

    input:
    tuple val(meta), path(mzml), path(idparquet)

    output:
    tuple val(meta),
        path(mzml),
        path("${meta.id}.percolator.idparquet"),
        emit: runs

    script:
    """
    PercolatorAdapter \
        -in "${idparquet}" \
        -out "${meta.id}.percolator.idparquet" \
        -use_subprocess true \
        -percolator_executable percolator \
        -score_type pep \
        -score:fdr ${params.psm_fdr} \
        -threads ${task.cpus}
    """
}
