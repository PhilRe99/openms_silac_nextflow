include { SDRF_PARSING } from '../../../modules/local/sdrf_parsing/main.nf'

workflow CREATE_INPUT_CHANNEL {

    take:
    ch_sdrf

    main:

    SDRF_PARSING(ch_sdrf)

    ch_openms = SDRF_PARSING.out.ch_sdrf_config_file

    ch_openms
        .splitCsv(header: true, sep: '\t')
        .map { row ->

            def mzml = file(
                "${params.spectra_dir ?: params.mzml_dir}/${row.Filename}",
                checkIfExists: true
            )

            def extension = mzml.extension.toLowerCase()

            if (!(extension in ['raw', 'mzml'])) {
                error "Unsupported input format: ${mzml.name}"
            }

            if (extension == 'raw' && !params.convert_raw) {
                error "RAW input ${mzml.name} requires --convert_raw"
            }

            def meta = [
                id: row.Filename.replaceFirst(/(?i)\.(raw|mzml)$/, ''),
                acquisition_method: row['Proteomics Data Acquisition Method'],
                label_type: row.Label,
                enzyme: row.Enzyme,
                fixed_modifications: row.FixedModifications,
                variable_modifications: row.VariableModifications,
                precursor_mass_tolerance: row.PrecursorMassTolerance,
                precursor_mass_tolerance_unit: row.PrecursorMassToleranceUnit,
                fragment_mass_tolerance: row.FragmentMassTolerance,
                fragment_mass_tolerance_unit: row.FragmentMassToleranceUnit,
                dissociation_method: row.DissociationMethod
            ]

            [meta, mzml]
        }
        .set { ch_runs }

    emit:
    runs = ch_runs
    openms = SDRF_PARSING.out.ch_sdrf_config_file
    experimental_design = SDRF_PARSING.out.ch_expdesign
}