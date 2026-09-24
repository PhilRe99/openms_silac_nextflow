include { SDRF_PARSING }       from '../modules/local/sdrf_parsing/main'

workflow CREATE_INPUT_CHANNEL {

    take:
    ch_sdrf

    main:

    SDRF_PARSING(ch_sdrf)

    ch_openms = SDRF_PARSING.out.ch_sdrf_config_file

    ch_openms
        .splitCsv(header: true, sep: '\t')
        .map { row ->

            def mzml = file("${params.mzml_dir}/${row.Filename}")

            def meta = [
                id: row.Filename.replaceFirst(/\.mzML$/, '')
            ]

            [meta, mzml]
        }
        .set { ch_runs }

    emit:
    runs = ch_runs
}