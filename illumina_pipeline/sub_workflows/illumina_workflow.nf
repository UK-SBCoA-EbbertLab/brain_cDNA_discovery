// Import Modules
include {DECOMPRESS ; TRIM_GALORE} from '../modules/trim_galore'
include {MAKE_STAR_INDEX ; STAR_MAPPING} from '../modules/STAR'
include {RSEQC} from '../modules/rseqc'
include {MULTIQC_GRCh38 ; MULTIQC_CHM13} from '../modules/multiqc'
include {SALMON_ALIGNMENT_MODE} from '../modules/salmon'


workflow ILLUMINA {

    take:
        ref
        gtf
        transcriptome
        housekeeping
        illumina_data
        overhang

    main:
        MAKE_STAR_INDEX(ref, gtf, overhang)

        DECOMPRESS(illumina_data)

        TRIM_GALORE(DECOMPRESS.out.id, DECOMPRESS.out.file_R1, DECOMPRESS.out.file_R2)
        
        STAR_MAPPING(TRIM_GALORE.out.id, TRIM_GALORE.out.trim_1, TRIM_GALORE.out.trim_2, MAKE_STAR_INDEX.out)

        SALMON_ALIGNMENT_MODE(STAR_MAPPING.out.id, STAR_MAPPING.out.bam_trans, transcriptome)

}

