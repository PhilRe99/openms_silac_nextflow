process PEAK_PICKER {

    tag "${meta.id}"

    input:
    tuple val(meta), path(mzml), path(config)

    output:
    tuple val(meta), path("prepared/${mzml.name}"), path(config), emit: runs

    script:
    """
    mkdir -p prepared

    PeakPickerHiRes \
        -in "${mzml}" \
        -out "prepared/${mzml.name}" \
        -threads ${task.cpus}
    """
}