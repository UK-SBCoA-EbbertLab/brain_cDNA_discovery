// Import Modules
include {DECOMPRESS ; TRIM_GALORE} from '../modules/trim_galore'
include {SALMON_ALIGNMENT_MODE} from '../modules/salmon'
include {MAKE_STAR_INDEX ; STAR_MAPPING_TRANSCRIPTOME} from '../modules/STAR'
include {DOWNSAMPLING_ILLUMINA} from '../modules/downsampling'
include {STRINGTIE_ILLUMINA} from '../modules/stringtie'


// TODO: Add RSEQC and FastQC to the quality control step for Illumina


workflow ILLUMINA {

    take:
        ref
        bambu_gtf
        transcriptome
        downsizes
        illumina_data

    main:
        MAKE_STAR_INDEX(ref, bambu_gtf)

        DECOMPRESS(illumina_data)

        TRIM_GALORE(DECOMPRESS.out.id, DECOMPRESS.out.file_R1, DECOMPRESS.out.file_R2)

        STAR_MAPPING_TRANSCRIPTOME(TRIM_GALORE.out.id, TRIM_GALORE.out.trim_1, TRIM_GALORE.out.trim_2, MAKE_STAR_INDEX.out)

        DOWNSAMPLING_ILLUMINA(STAR_MAPPING_TRANSCRIPTOME.out.id, STAR_MAPPING_TRANSCRIPTOME.out.bam_trans, STAR_MAPPING_TRANSCRIPTOME.out.bam_gen, downsizes)

        SALMON_ALIGNMENT_MODE(DOWNSAMPLING_ILLUMINA.out.id, DOWNSAMPLING_ILLUMINA.out.downsample, DOWNSAMPLING_ILLUMINA.out.bam_salmon, transcriptome)

        STRINGTIE_ILLUMINA(DOWNSAMPLING_ILLUMINA.out.id, DOWNSAMPLING_ILLUMINA.out.downsample, DOWNSAMPLING_ILLUMINA.out.bam_stringtie, DOWNSAMPLING_ILLUMINA.out.bai_stringtie, bambu_gtf)


}
