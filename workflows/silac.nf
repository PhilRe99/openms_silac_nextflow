include { CREATE_INPUT_CHANNEL } from '../subworkflows/local/create_input_channel/main'
include { RESOLVE_SILAC_CONFIG } from '../modules/local/resolve_silac_config/main'
include { OPENMS_PEAK_PICKER } from '../modules/local/openms/openms_peak_picker/main'
include { GENERATE_DECOY_DATABASE } from '../modules/local/openms/generate_decoy_database/main'
include { COMET } from '../modules/local/openms/comet/main'
include { PEPTIDE_INDEXER } from '../modules/local/openms/peptide_indexer/main'
include { PSM_FEATURE_EXTRACTOR } from '../modules/local/openms/psm_feature_extractor/main'
include { PERCOLATOR } from '../modules/local/openms/percolator/main'
include { MS1_LABELED_WORKFLOW } from '../modules/local/openms/ms1_labeled_workflow/main'

workflow SILAC {

    main:

    CREATE_INPUT_CHANNEL(
        Channel.value(file(params.input))
    )

    //extracting silac specific info from metadatat into config file
    RESOLVE_SILAC_CONFIG(
        CREATE_INPUT_CHANNEL.out.runs,
        CREATE_INPUT_CHANNEL.out.openms,
        CREATE_INPUT_CHANNEL.out.experimental_design,
        Channel.value(file(params.silac_chemistry, checkIfExists: true))
    )

    //expanding meta for the needed additional info
    RESOLVE_SILAC_CONFIG.out.runs.map { meta, mzml, config ->
        def cfg = new groovy.json.JsonSlurperClassic().parse(config.toFile())
        def final_meta = meta + [
            variable_modifications: cfg.variable_modifications,
            binary_modifications: cfg.binary_modifications,
            silac_labels: cfg.ffm_labels
        ]
        tuple(final_meta, mzml)
    }.set { resolved_runs }

    //only accepts one shared chemistry for now (no mixed plex)
    resolved_runs
    .collect(flat: false)
    .map { runs ->
        def chemistries = runs.collect { run -> run[0].silac_labels }.unique()

        assert chemistries.size() == 1 :
            "Expected one SILAC chemistry, found ${chemistries}"

        runs
    }
    .flatMap { runs -> runs }
    .set { validated_runs }

    OPENMS_PEAK_PICKER(
        validated_runs
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
        GENERATE_DECOY_DATABASE.out.database
    )

    PSM_FEATURE_EXTRACTOR(
        PEPTIDE_INDEXER.out.runs
    )

    PERCOLATOR(
        PSM_FEATURE_EXTRACTOR.out.runs
    )

    //pooling for MS1labeledworkflow
    PERCOLATOR.out.runs.collect(flat: false).map { runs ->
        def chemistries = runs.collect { run -> run[0].silac_labels }.unique()
        tuple(chemistries[0], runs)
    }.set { ms1_runs }
    

    //restructuring collected runs into form the MS1 workflow needs
    //run[1] = mzML
    //run[2] = idXML
    ms1_runs.map { labels, runs ->
        tuple(
            labels,
            runs.collect { run -> run[1] },
            runs.collect { run -> run[2] }
        )
    }.set { ms1_input }

    MS1_LABELED_WORKFLOW(
    ms1_input,
    CREATE_INPUT_CHANNEL.out.experimental_design,
    GENERATE_DECOY_DATABASE.out.database
    )
}
