include { CREATE_INPUT_CHANNEL } from '../subworkflows/local/create_input_channel/main'
include {RESOLVE_SILAC_CONFIG} from '../modules/local/resolve_silac_config/main'
include { PEAK_PICKER } from '../modules/local/peak_picker/main'

workflow SILAC {

    main:

    CREATE_INPUT_CHANNEL(
        Channel.value(file(params.input))
    )

    RESOLVE_SILAC_CONFIG(
        CREATE_INPUT_CHANNEL.out.runs,
        CREATE_INPUT_CHANNEL.out.openms,
        CREATE_INPUT_CHANNEL.out.experimental_design
    )

    PEAK_PICKER(
    RESOLVE_SILAC_CONFIG.out.runs
    )

    PEAK_PICKER.out.runs.view()

}