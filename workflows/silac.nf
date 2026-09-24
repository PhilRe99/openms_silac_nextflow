include { CREATE_INPUT_CHANNEL } from '../subworkflows/local/create_input_channel/main'

workflow SILAC {


    main:

    CREATE_INPUT_CHANNEL(
        Channel.value(file(params.input))
    )

    CREATE_INPUT_CHANNEL.out.runs.view()

}