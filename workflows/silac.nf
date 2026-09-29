include { CREATE_INPUT_CHANNEL } from '../subworkflows/local/create_input_channel/main'
include { RESOLVE_SILAC_CONFIG } from '../modules/local/resolve_silac_config/main'
include { OPENMS_PEAK_PICKER } from '../modules/local/openms/openms_peak_picker/main'
include { GENERATE_DECOY_DATABASE } from '../modules/local/openms/generate_decoy_database/main'
include { COMET } from '../modules/local/openms/comet/main'
include { PEPTIDE_INDEXER } from '../modules/local/openms/peptide_indexer/main'
include { PSM_FEATURE_EXTRACTOR } from '../modules/local/openms/psm_feature_extractor/main'
include { PERCOLATOR } from '../modules/local/openms/percolator/main'

workflow SILAC {

    main:

    CREATE_INPUT_CHANNEL(
        Channel.value(file(params.input))
    )

    RESOLVE_SILAC_CONFIG(
        CREATE_INPUT_CHANNEL.out.runs,
        CREATE_INPUT_CHANNEL.out.openms.first(),
        CREATE_INPUT_CHANNEL.out.experimental_design.first()
    )

    RESOLVE_SILAC_CONFIG.out.runs
        .map { meta, mzml, config ->
            def cfg = new groovy.json.JsonSlurperClassic().parse(config.toFile())
            def resolved_meta = meta + [
                labels: cfg.labels,
                label_channels: cfg.label_channels,
                label_modifications: cfg.label_modifications,
                variable_modifications: cfg.variable_modifications,
                binary_modifications: cfg.binary_modifications,
                silac_labels: cfg.ffm_labels
            ]
            tuple(final_meta, mzml)
        }
        .set { resolved_runs }

    OPENMS_PEAK_PICKER(
        resolved_runs
    )

    GENERATE_DECOY_DATABASE(
        Channel.value(file(params.database, checkIfExists: true))
    )

    OPENMS_PEAK_PICKER.out.runs
        .combine(GENERATE_DECOY_DATABASE.out.database)
        .set { comet_input }

    COMET(comet_input)

    PEPTIDE_INDEXER(
        COMET.out.runs,
        GENERATE_DECOY_DATABASE.out.database.first()
    )

    PSM_FEATURE_EXTRACTOR(
        PEPTIDE_INDEXER.out.runs
    )

    PERCOLATOR(
        PSM_FEATURE_EXTRACTOR.out.runs
    )

    PERCOLATOR.out.runs.view()

}
