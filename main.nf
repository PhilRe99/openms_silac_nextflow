nextflow.enable.dsl = 2

include {SILAC} from './workflows/silac'

workflow {
    SILAC()
}