// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2


// Pipeline parameter default values, can be modified by user when calling pipeline on command line (e.g. --data_fq sample_1.fastq) ##
params.reads_fq = '/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/*.fastq'
params.reads_txt = '/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/*.txt'
params.ref = '../ont_references/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.gtf = '../ont_references/Homo_sapiens.GRCh38.104.gtf'
params.bed = "../ont_references/hg38_ncbi_08_29_2021.bed"
params.housekeep = "../ont_references/hg38.HouseKeepingGenes.bed"


log.info """\
 RNA-SEQ ONT CORRECTION AND DISCOVERY PIPELINE
 ===============================================
 fastq files                   : ${params.reads_fq}
 sequencing summary files      : ${params.reads_txt}
 reference genome              : ${params.ref}
 reference annotation          : ${params.gtf}
 reference bed                 : ${params.bed}
 reference housekeeping bed    : ${params.housekeep}
 """


// Import Modules
include {MAKE_FAI} from '../ont_modules/make_fai'
include {MAKE_INDEX} from '../ont_modules/make_index'
include {PYCHOPPER} from '../ont_modules/pychopper'
include {MAPPING ; MAPPING_QC ; MAPPING_QC_DIRECT_RNA} from '../ont_modules/mapping'
include {PYCOQC} from '../ont_modules/pycoqc'
include {RSEQC} from '../ont_modules/rseqc'
include {DOWNSAMPLING} from '../ont_modules/downsampling'
include {BAMBU} from '../ont_modules/bambu'

// Define initial files and channels
reads_fq = Channel.fromPath(params.reads_fq).map { file -> tuple(file.baseName, file) }
reads_txt = Channel.fromPath(file(params.reads_txt))
ref = file(params.ref)
gtf = file(params.gtf)
bed = file(params.bed)
housekeep = file(params.housekeep)
// downsizes = Channel.from(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
downsizes = Channel.from(1.00)


workflow {

    MAKE_FAI(ref)

    MAKE_INDEX(ref)

    // MAPPING_QC(reads_fq, MAKE_INDEX.out)

    // PYCOQC(MAPPING_QC.out.id, reads_txt, MAPPING_QC.out.bam, MAPPING_QC.out.bai)

    PYCHOPPER(reads_fq)

    MAPPING(PYCHOPPER.out.id, PYCHOPPER.out.fastq, MAKE_INDEX.out)

   //  RSEQC(MAPPING.out.id, MAPPING.out.bam, MAPPING.out.bai, housekeep, bed)

    DOWNSAMPLING(MAPPING.out.id, MAPPING.out.bam, downsizes)

    BAMBU(DOWNSAMPLING.out.id, DOWNSAMPLING.out.bam.groupTuple(), DOWNSAMPLING.out.bai.groupTuple(), ref, gtf, MAKE_FAI.out)
}
