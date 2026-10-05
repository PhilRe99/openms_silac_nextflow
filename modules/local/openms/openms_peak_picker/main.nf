process OPENMS_PEAK_PICKER {

    tag "${meta.id}"

    container params.openms_container

    input:
    tuple val(meta), path(mzml)

    output:
    tuple val(meta), path("prepared/${mzml.name}"), emit: runs

    script:
    """
    mkdir -p prepared

    PeakPickerHiRes \
        -in "${mzml}" \
        -out "prepared/${mzml.name}" \
        -threads ${task.cpus}
    """
}
