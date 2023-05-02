// Import Modules
include {MAKE_FAI} from '../modules/make_fai'
include {MAKE_INDEX_dRNA} from '../modules/make_index'
include {PYCHOPPER} from '../modules/pychopper'
include {MINIMAP2_dRNA} from '../modules/minimap2'
include {DOWNSAMPLING_NANOPORE} from '../modules/downsampling'
include {BAMBU_PREP_DOWN; BAMBU_SPIKE} from '../modules/bambu'
include {STRINGTIE_ONT_dRNA} from '../modules/stringtie'
include {MAKE_TRANSCRIPTOME} from '../modules/make_transcriptome'


workflow SPIKE_IN_dRNA {

    take:
        ref
        annotation
        downsizes
        ont_reads_txt
        ont_reads_fq

    main:
        
        downsizes = Channel.from(1.00)

        MAKE_FAI(ref)
        MAKE_INDEX_dRNA(ref)
        MINIMAP2_dRNA(ont_reads_fq, MAKE_INDEX_dRNA.out)
        DOWNSAMPLING_NANOPORE(MINIMAP2_dRNA.out.id, MINIMAP2_dRNA.out.bam_mapped, 1.00)
        BAMBU_PREP_DOWN(DOWNSAMPLING_NANOPORE.out.bam, DOWNSAMPLING_NANOPORE.out.bai, DOWNSAMPLING_NANOPORE.out.downsample, ref, annotation, MAKE_FAI.out)
        BAMBU_SPIKE(BAMBU_PREP_DOWN.out.groupTuple(), ref, annotation, MAKE_FAI.out)
        STRINGTIE_ONT_dRNA(MINIMAP2_dRNA.out.id, MINIMAP2_dRNA.out.bam_mapped, MINIMAP2_dRNA.out.bai_mapped, annotation)
        MAKE_TRANSCRIPTOME(ref, MAKE_FAI.out, annotation)


}
