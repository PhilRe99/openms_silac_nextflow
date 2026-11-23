process THERMORAWFILEPARSER {

    tag "${meta.id}"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("${meta.id}.mzML"), emit: runs

    script:
    """
    ThermoRawFileParser.sh \\
        -i="${raw}" \\
        -b="${meta.id}.mzML" \\
        -f=2
    """
}
