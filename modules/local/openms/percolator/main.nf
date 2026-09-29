process PERCOLATOR {

    tag "${meta.id}"

    input:
    tuple val(meta), path(mzml), path(idxml), path(config)

    output:
    tuple val(meta),
        path(mzml),
        path("${meta.id}.percolator.idXML"),
        path(config),
        emit: runs

    script:
    """
    PercolatorAdapter \
        -in "${idxml}" \
        -out "${meta.id}.percolator.idXML" \
        -use_subprocess true \
        -percolator_executable /usr/share/OpenMS/THIRDPARTY/Percolator/percolator \
        -score_type pep \
        -score:fdr ${params.psm_fdr}
        -threads ${task.cpus}
    """
}