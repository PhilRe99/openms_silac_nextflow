process THERMORAWFILEPARSER {
    tag "${meta.id}"

    container 'biocontainers/thermorawfileparser:1.4.5--h05cac1d_1'

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("${meta.id}.mzML"), emit: runs

    script:
    """
    ThermoRawFileParser.sh \\
        -i="${raw}" \\
        -b="${meta.id}.mzML" \\
        -f=2 \\
        -p
    """
}