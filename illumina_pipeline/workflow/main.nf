// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2


log.info """
 cDNA SEQUENCING PIPELINE FOR ILLUMINA
 ===========================================================================
 illumina fastq files                           : ${params.illumina_data}
 reference genome                               : ${params.ref}
 reference annotation                           : ${params.annotation}
 reference transcriptome                        : ${params.transcriptome}
 housekeeping genes 3' bias assessment          : ${params.housekeeping}
 overhang (combined read size - 1)              : ${params.overhang}
"""


// Import Workflows
include {ILLUMINA} from '../sub_workflows/illumina_workflow'

// Define initial files and channels
illumina_data = Channel.fromFilePairs(params.illumina_data, flat:true)
ref = file(params.ref)
housekeeping = file(params.housekeeping)
annotation = file(params.annotation)
overhang = Channel.value(params.overhang)
transcriptome = file(params.transcriptome)

workflow {

    ILLUMINA(ref, annotation, transcriptome, housekeeping, illumina_data, overhang)
}
