process SDRF_PARSING {
     
    tag "${sdrf.simpleName}"
    
    container 'openms-silac/sdrf-pipelines:dea8ba9-procps'

    input:
    path sdrf

    output:
    path "metadata/${sdrf.name}", emit: source_sdrf
    path "${sdrf.baseName}_openms_design.tsv", emit: ch_expdesign
    path "${sdrf.baseName}_config.tsv"       , emit: ch_sdrf_config_file

    script:
    def extensionArgs = params.convert_raw.toString().toBoolean() ? '--extension_convert raw:mzML' : ''

    """
    parse_sdrf validate-sdrf   --sdrf_file "$sdrf"   --template ms-proteomics

    mkdir -p metadata
    cp "${sdrf}" "metadata/${sdrf.name}"

    parse_sdrf convert-openms -s "$sdrf" ${extensionArgs}

    mv openms.tsv ${sdrf.baseName}_config.tsv
    mv experimental_design.tsv ${sdrf.baseName}_openms_design.tsv
    """
}
