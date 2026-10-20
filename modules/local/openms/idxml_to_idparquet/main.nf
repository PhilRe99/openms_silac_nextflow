process IDXML_TO_IDPARQUET {

    tag "${meta.id}"

    container params.openms_container

    input:
    tuple val(meta), path(mzml), path(identifications)

    output:
    tuple val(meta), path(mzml),
        path("${meta.id}.percolator.idparquet", type: 'dir'), emit: runs

    script:
    """
    IDFileConverter \
        -in "${identifications}" \
        -out "${meta.id}.percolator.idparquet" \
        -threads ${task.cpus}
    """
}
