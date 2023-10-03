// Import Modules
include {FILTER_UNIQUE_BAM} from '../modules/filter_unique_bam'
include {SALMON_ALIGNMENT_MODE} from '../modules/salmon'


workflow FILTER_ILLUMINA_UNIQUE {

    take:
        bam
        transcriptome

    main:
        
        FILTER_UNIQUE_BAM(bam)
        
        SALMON_ALIGNMENT_MODE(FILTER_UNIQUE_BAM.out.id, FILTER_UNIQUE_BAM.out.bam, transcriptome)

}

