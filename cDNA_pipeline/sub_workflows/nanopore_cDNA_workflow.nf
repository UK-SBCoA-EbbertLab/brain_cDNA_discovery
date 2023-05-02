// Import Modules
include {MAKE_FAI} from '../modules/make_fai'
include {MAKE_INDEX_cDNA} from '../modules/make_index'
include {CHM13_GTF; CHM13_GTF_ERCC} from '../modules/chm13_gff3_to_gtf'
include {PYCHOPPER} from '../modules/pychopper'
include {PYCOQC} from '../modules/pycoqc'
include {MINIMAP2_cDNA} from '../modules/minimap2'
include {RSEQC} from '../modules/rseqc'
include {BAMBU_PREP_DISCOVERY; BAMBU_PREP_QUANT; BAMBU_DISCOVERY; BAMBU_QUANT} from '../modules/bambu'
include {STRINGTIE_ONT_cDNA} from '../modules/stringtie'
include {GFFCOMPARE} from '../modules/gffcompare'
include {MAKE_TRANSCRIPTOME} from '../modules/make_transcriptome'
include {MULTIQC_GRCh38 ; MULTIQC_CHM13} from '../modules/multiqc'


workflow NANOPORE_cDNA {

    take:
        ref
        annotation
        housekeeping
        ont_reads_txt
        ont_reads_fq
        ercc
        cdna_kit
        multiqc_config
        NDR
        track_reads

    main:
        MAKE_FAI(ref)
        MAKE_INDEX_cDNA(ref)
        PYCHOPPER(ont_reads_fq, ont_reads_txt, cdna_kit)
        MINIMAP2_cDNA(PYCHOPPER.out.id, PYCHOPPER.out.fastq,  MAKE_INDEX_cDNA.out, PYCHOPPER.out.txt)
        PYCOQC(MINIMAP2_cDNA.out.id, MINIMAP2_cDNA.out.fastq, MINIMAP2_cDNA.out.txt, MINIMAP2_cDNA.out.bam_all, MINIMAP2_cDNA.out.bai_all)

        if (params.is_chm13 == true)
        {
            MULTIQC_CHM13(MINIMAP2_cDNA.out.QC_out.collect(), PYCOQC.out.multiQC.collect(), PYCHOPPER.out.multiQC.collect(), multiqc_config)

            if (params.ercc == "None") 
            { 
                CHM13_GTF(annotation)
                annotation = CHM13_GTF.out.collect()
            }
            
            else 
            {
                CHM13_GTF_ERCC(annotation, ercc)
                annotation = CHM13_GTF_ERCC.out.collect()
            }
        }

        else
        {
            RSEQC(MINIMAP2_cDNA.out.bam_mapped.collect(), MINIMAP2_cDNA.out.bai_mapped.collect(), housekeeping)
            MULTIQC_GRCh38(MINIMAP2_cDNA.out.QC_out.collect(), PYCOQC.out.multiQC.collect(), PYCHOPPER.out.multiQC.collect(), RSEQC.out.multiQC.collect(), multiqc_config)
        }
        
        if (params.is_discovery == true)
        {

            BAMBU_PREP_DISCOVERY(MINIMAP2_cDNA.out.bam_mapped, MINIMAP2_cDNA.out.bai_mapped, ref, annotation, MAKE_FAI.out, NDR, track_reads)
            BAMBU_DISCOVERY(BAMBU_PREP_DISCOVERY.out.collect(), ref, annotation, MAKE_FAI.out, NDR, track_reads)
            new_annotation = BAMBU_DISCOVERY.out.gtf
            GFFCOMPARE(new_annotation, annotation)

        }

        else
        {
            
            BAMBU_PREP_QUANT(MINIMAP2_cDNA.out.bam_mapped, MINIMAP2_cDNA.out.bai_mapped, ref, annotation, MAKE_FAI.out)
            BAMBU_QUANT(BAMBU_PREP_QUANT.out.collect(), ref, annotation, MAKE_FAI.out)
            new_annotation = BAMBU_QUANT.out.gtf

        }

        STRINGTIE_ONT_cDNA(MINIMAP2_cDNA.out.id, MINIMAP2_cDNA.out.bam_mapped, MINIMAP2_cDNA.out.bai_mapped, new_annotation)
        MAKE_TRANSCRIPTOME(ref, MAKE_FAI.out, new_annotation)
        
}
