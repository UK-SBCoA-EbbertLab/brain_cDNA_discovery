// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2


log.info """
 RNA-SEQ ONT CORRECTION AND DISCOVERY PIPELINE
 ===============================================
 nanopore fastq files                   : ${params.ont_reads_fq}
 nanopore sequencing summary files      : ${params.ont_reads_txt}
 illumina fastq files                   : ${params.ill_reads_fq}
 reference genome                       : ${params.ref}
 reference annotation                   : ${params.gtf}
 """


// Import Workflows
include {ILLUMINA} from '../sub_workflows/illumina_workflow'
include {NANOPORE_CDNA} from '../sub_workflows/nanopore_workflow'


// Define initial files and channels
ont_reads_fq = Channel.fromPath(params.ont_reads_fq).map { file -> tuple(file.baseName, file) }
ont_reads_txt = Channel.fromPath(file(params.ont_reads_txt))
ill_reads_fq = Channel.fromFilePairs(params.ill_reads_fq, flat:true)
ref = file(params.ref)
gtf = file(params.gtf)
bed = file(params.bed)
housekeep = file(params.housekeep)
downsizes = Channel.from(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
illumina_data = Channel.fromFilePairs(params.illumina_data, flat:true)


// Make sure ONT sequencing summary and fastq files are in the same order
ont_reads_txt = ont_reads_txt.toSortedList( { a, b -> a.baseName <=> b.baseName } ).flatten()
ont_reads_fq = ont_reads_fq.toSortedList( { a, b -> a[0] <=> b[0] } ).flatten().buffer(size:2)

workflow {
    
    NANOPORE_CDNA(ref, gtf, bed, housekeep, downsizes, ont_reads_txt, ont_reads_fq)
    ILLUMINA(ref, NANOPORE_CDNA.out[0], NANOPORE_CDNA.out[1], downsizes, illumina_data)

}

