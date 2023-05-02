// Import Modules
include {MAKE_FAI} from '../modules/make_fai'
include {MAKE_INDEX_dRNA} from '../modules/make_index'
include {CHM13_GTF; FIX_EXTENDED_ANNOTATION} from '../modules/chm13_gff3_to_gtf'
include {MINIMAP2_dRNA} from '../modules/minimap2'
include {RSEQC} from '../modules/rseqc'
include {BAMBU_dRNA} from '../modules/bambu'
include {STRINGTIE_ONT_dRNA} from '../modules/stringtie'
include {GFFCOMPARE} from '../modules/gffcompare'


// TODO: MAKE QUALITY CONTROL WORK (RSEQC + PYCOQC + MULTIQC) ONCE WE HAVE FINAL DATA

workflow NANOPORE_dRNA {

    take:
        ref
        annotation
        downsizes
        ont_reads_txt
        ont_reads_fq
    
    main:

        MAKE_FAI(ref)

        MAKE_INDEX_dRNA(ref)

        MINIMAP2_dRNA(ont_reads_fq, MAKE_INDEX_dRNA.out)

        if (params.is_chm13 == true)
        {
            CHM13_GTF(annotation)
            annotation = CHM13_GTF.out
        }


        RSEQC(MINIMAP2_dRNA.out.id, MINIMAP2_dRNA.out.bam_all, MINIMAP2_dRNA.out.bai_all)

        BAMBU_dRNA(MINIMAP2_dRNA.out.bam_mapped.collect(), MINIMAP2_dRNA.out.bai_mapped.collect(), ref, annotation, MAKE_FAI.out)

        STRINGTIE_ONT_dRNA(MINIMAP2_dRNA.out.id, MINIMAP2_dRNA.out.bam_mapped, MINIMAP2_dRNA.out.bai_mapped, BAMBU_dRNA.out.gtf)
    
        GFFCOMPARE(BAMBU_dRNA.out.gtf, annotation)
}
