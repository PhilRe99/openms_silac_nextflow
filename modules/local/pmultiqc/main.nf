process PMULTIQC {

    input:
    path 'results/*'

    output:
    path '*.html', emit: ch_pmultiqc_report
    path '*.db', optional: true, emit: ch_pmultiqc_db
    path 'versions.yml', emit: versions
    path '*_data', emit: data

    script:
    def args = task.ext.args ?: ''
    """
    multiqc -f ${args} ./results -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pmultiqc: \$(multiqc --pmultiqc-version | sed -e "s/pmultiqc, version //g")
    END_VERSIONS
    """
}
