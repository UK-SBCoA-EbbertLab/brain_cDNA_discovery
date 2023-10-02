// Import Modules
include {FILTER_UNIQUE_BAM} from '../modules/filter_bam_unique'
include {SALMON_ALIGNMENT_MODE} from '../modules/salmon'


workflow ILLUMINA {

    take:
        transcriptome
        bam

    main:
        
        FILTER_UNIQUE_BAM(bam)
        
        SALMON_ALIGNMENT_MODE(FILTER_UNIQUE_BAM.out.id, FILTER_UNIQUE_BAM.out.bam_trans, transcriptome)

}

