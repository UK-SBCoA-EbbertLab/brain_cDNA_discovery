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
include {FILTER_ILLUMINA_UNIQUE} from '../sub_workflows/illumina_workflow_from_BAM.nf'

// Define initial files and channels

if (params.bam == "None") {
    
    illumina_data = Channel.fromFilePairs(params.illumina_data, flat:true)
    ref = file(params.ref)
    housekeeping = file(params.housekeeping)
    annotation = file(params.annotation)
    overhang = Channel.value(params.overhang)

    }
else {

    bam = Channel.fromPath(params.bam).map { file -> tuple(file.baseName, file) }

}

transcriptome = file(params.transcriptome)

workflow {
    
    if (params.bam == "None") {
    
    ILLUMINA(ref, annotation, transcriptome, housekeeping, illumina_data, overhang)

    }

    else {
    
    FILTER_ILLUMINA_UNIQUE(bam, transcriptome)

    }
}
