// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2


// Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --data_fq sample_1.fastq) ##
params.ref = '../illumina_references/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.gtf = '../illumina_references/Homo_sapiens.GRCh38.105.gtf'
params.data_dir = '/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_illumina_data/'
params.fastq_pattern = '*_R{1,2}*.fastq.gz'
params.multiqc_config = "../illumina_references/multiqc_config.yaml"
params.transcriptome = '../illumina_references/Homo_sapiens.GRCh38.transcriptome.fa'


log.info """\
 ILLUMINA RNA-SEQ PIPELINE FOR CDNA COMPARISON PROJECT
 CREATED BY BERNARDO AGUZZOLI HEBERLE - EBBERT LAB
 ===============================================
 data directory with fastq files                   : ${params.data_dir}
 pattern of paired fastq files                     : ${params.fastq_pattern}
 reference genome                                  : ${params.ref}
 reference annotation                              : ${params.gtf}
 salmon transcriptome reference                    : ${params.transcriptome}
 multiqc report configuration                      : ${params.multiqc_config}
 """


// Import Modules
include {DECOMPRESS ; TRIM_GALORE} from '../illumina_modules/trim_galore'
include {SALMON_MAPPING_MODE ; SALMON_ALIGNMENT_MODE} from '../illumina_modules/salmon'
include {MAKE_STAR_INDEX ; STAR_MAPPING_TRANSCRIPTOME} from '../illumina_modules/STAR'
include {DOWNSAMPLING} from '../illumina_modules/downsampling'
include {MULTIQC} from '../illumina_modules/multiqc'


// Define initial files and channels
fastq_files = params.data_dir + params.fastq_pattern
downsizes = Channel.from(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
annotation = Channel.value(file(params.gtf))
reference = Channel.value(file(params.ref))
fastq_list = Channel.fromFilePairs(fastq_files, flat:true)
multiqc = Channel.value(file(params.multiqc_config))
transcriptome = Channel.value(file(params.transcriptome))

// Workflow Implementation

workflow {
    
    DECOMPRESS(fastq_list)

    TRIM_GALORE(DECOMPRESS.out.id, DECOMPRESS.out.file_R1, DECOMPRESS.out.file_R2)

    MAKE_STAR_INDEX(reference, annotation)

    STAR_MAPPING_TRANSCRIPTOME(TRIM_GALORE.out.id, TRIM_GALORE.out.trim_1, TRIM_GALORE.out.trim_2, MAKE_STAR_INDEX.out)

    DOWNSAMPLING(STAR_MAPPING_TRANSCRIPTOME.out.id, STAR_MAPPING_TRANSCRIPTOME.out.bam, downsizes)

    SALMON_ALIGNMENT_MODE(DOWNSAMPLING.out.id, DOWNSAMPLING.out.downsample, DOWNSAMPLING.out.bam, transcriptome)
        
    //MULTIQC(TRIM_GALORE.out.QC_1.collect(), TRIM_GALORE.out.QC_2.collect(), MAPPING.out.QC_3.collect(), MAPPING.out.QC_4.collect(), multiqc)

}
