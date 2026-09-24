process SDRF_PARSING {
     
    tag "${sdrf.simpleName}"
    
    container 'openms-silac/sdrf-pipelines:dea8ba9'

    input:
    path sdrf

    output:
    path "${sdrf.baseName}_openms_design.tsv", emit: ch_expdesign
    path "${sdrf.baseName}_config.tsv"       , emit: ch_sdrf_config_file

    script:
    """
    parse_sdrf validate-sdrf   --sdrf_file "$sdrf"   --template ms-proteomics
    parse_sdrf convert-openms -s "$sdrf"

    mv openms.tsv ${sdrf.baseName}_config.tsv
    mv experimental_design.tsv ${sdrf.baseName}_openms_design.tsv
    """
}