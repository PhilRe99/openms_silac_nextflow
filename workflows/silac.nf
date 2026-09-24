include { SDRF_PARSING }       from '../modules/local/sdrf_parsing/main'
include { CREATE_INPUT_CHANNEL } from '../subworkflows/local/create_input_channel/main'

workflow SILAC {


    main:

    SDRF_PARSING(
        file(params.input)
    )

    CREATE_INPUT_CHANNEL(
        SDRF_PARSING.out.ch_sdrf_config_file
    )

    CREATE_INPUT_CHANNEL.out.runs.view()
    
}