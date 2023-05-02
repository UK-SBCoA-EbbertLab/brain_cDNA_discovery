// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2


// Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --data_fq sample_1.fastq) ##
params.illumina_data = '/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_illumina_data/*_R{1,2}*.fastq.gz'
params.ont_reads_fq = '/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/*.fastq'
params.ont_reads_txt = '/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/*.txt'
params.ill_reads_fq = '/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_illumina_data/*_R{1,2}*.fastq.gz'
params.ref = '../ont_references/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.gtf = '../ont_references/Homo_sapiens.GRCh38.104.gtf'
params.bed = "../ont_references/hg38_ncbi_08_29_2021.bed"
params.housekeep = "../ont_references/hg38.HouseKeepingGenes.bed"

log.info """
 RNA-SEQ ONT CORRECTION AND DISCOVERY PIPELINE
 ===============================================
 nanopore fastq files                   : ${params.ont_reads_fq}
 nanopore sequencing summary files      : ${params.ont_reads_txt}
 illumina fastq files                   : ${params.ill_reads_fq}
 reference genome                       : ${params.ref}
 reference annotation                   : ${params.gtf}
 reference bed                          : ${params.bed}
 reference housekeeping bed             : ${params.housekeep}
 """


// Import Modules
include {MAKE_FAI} from '../integrated_modules/make_fai'
include {MAKE_INDEX} from '../integrated_modules/make_index'
include {PYCHOPPER} from '../integrated_modules/pychopper'
include {SEQ_SUMMARY} from '../integrated_modules/fix_sequencing_summary'
include {MINIMAP2_CDNA} from '../integrated_modules/minimap2'
include {PYCOQC} from '../integrated_modules/pycoqc'
// include {RSEQC} from '../integrated_modules/rseqc'
include {DOWNSAMPLING_NANOPORE} from '../integrated_modules/downsampling'
include {BAMBU_DISCOVERY; BAMBU_QUANT} from '../integrated_modules/bambu'
// include {STRINGTIE_QUANT} from '../integrated_modules/stringtie'
include {WEBSITE_ANNOTATIONS} from '../integrated_modules/expression_and_annotation_editor'
include {MAKE_TRANSCRIPTOME} from '../integrated_modules/make_transcriptome'
include {DECOMPRESS ; TRIM_GALORE} from '../illumina_modules/trim_galore'
include {SALMON_MAPPING_MODE ; SALMON_ALIGNMENT_MODE} from '../illumina_modules/salmon'
include {MAKE_STAR_INDEX ; STAR_MAPPING_TRANSCRIPTOME} from '../illumina_modules/STAR'
include {DOWNSAMPLING} from '../illumina_modules/downsampling'
include {MULTIQC} from '../illumina_modules/multiqc'



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


workflow NANOPORE_CDNA {

    take:
        ref
        gtf
        bed
        housekeep
        downsizes
        ont_reads_txt
        ont_reads_fq

    main:
        MAKE_FAI(ref)
        MAKE_INDEX(ref)
        PYCHOPPER(ont_reads_fq, ont_reads_txt)
        SEQ_SUMMARY(PYCHOPPER.out.id, PYCHOPPER.out.fastq, PYCHOPPER.out.txt)
        MINIMAP2_CDNA(SEQ_SUMMARY.out.id, SEQ_SUMMARY.out.fastq, SEQ_SUMMARY.out.txt, MAKE_INDEX.out)
        // PYCOQC(MINIMAP2_CDNA.out.txt, MINIMAP2_CDNA.out.bam_all, MINIMAP2_CDNA.out.bai_all)
        // RSEQC(MINIMAP2_CDNA.out.id, MINIMAP2_CDNA.out.bam_mapped, MINIMAP2_CDNA.out.bai_mapped, housekeep, bed)
        BAMBU_DISCOVERY(MINIMAP2_CDNA.out.bam_mapped.collect(), MINIMAP2_CDNA.out.bai_mapped.collect(), ref, gtf, MAKE_FAI.out)
        WEBSITE_ANNOTATIONS(gtf, BAMBU_DISCOVERY.out.gtf, BAMBU_DISCOVERY.out.counts)
        DOWNSAMPLING_NANOPORE(MINIMAP2_CDNA.out.id, MINIMAP2_CDNA.out.bam_mapped, downsizes)
        BAMBU_QUANT(DOWNSAMPLING_NANOPORE.out.bam.groupTuple(), DOWNSAMPLING_NANOPORE.out.bai.groupTuple(), ref, BAMBU_DISCOVERY.out.gtf, MAKE_FAI.out)
        // STRINGTIE_QUANT(DOWNSAMPLING.out.id, DOWNSAMPLING.out.bam.groupTuple(), DOWNSAMPLING.out.bai.groupTuple(), ref, BAMBU_DISCOVERY.out.gtf, MAKE_FAI.out))
        MAKE_TRANSCRIPTOME(ref, MAKE_FAI.out, BAMBU_DISCOVERY.out.gtf)

    emit:
        MAKE_TRANSCRIPTOME.out
        
}


workflow ILLUMINA {

    take:
        ref
        gtf_bambu
        downsizes
        illumina_reads
        

    
    main:
        NANOPORE_CDNA.out.view()
        downsizes.view()
        illumina_reads.view()

        MAKE_STAR_INDEX(ref, gtf_bambu)

}

workflow {
    
    take:
        ref
        gtf
        bed
        housekeep
        downsizes
        ont_reads_txt
        ont_reads_fq
        illumina_reads


    main:
        NANOPORE_CDNA(ref, gtf, bed, housekeep, downsizes, ont_reads_txt, ont_reads_fq)
        ILLUMINA(ref, NANOPORE_CDNA.out, downsizes, illumina_reads)

}


