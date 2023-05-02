params.data_fq = 'raw_data/*.fastq'
params.data_txt = 'raw_data/*.txt'
params.ref = 'references/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.gtf = 'references/Homo_sapiens.GRCh38.104.gtf'
params.name = "sample"
params.bed = "references/hg38_ncbi_08_29_2021.bed"
params.housekeep = "references/hg38.HouseKeepingGenes.bed"

data_fq = file(params.data_fq)
data_txt = file(params.data_txt)
ref = file(params.ref)
gtf = file(params.gtf)
name = params.name
bed = file(params.bed)
housekeep = file(params.housekeep)


order_1 = Channel.from(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
order_2 = Channel.from(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
downsizes = Channel.from(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
downsizes_2 = Channel.from(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
annotation = Channel.value(gtf)
reference = Channel.value(ref) 
raw_fq = Channel.value(data_fq)
raw_txt = Channel.value(data_txt)
bed_file = Channel.from(bed)
housekeep_file = Channel.from(housekeep)

process make_fai {

    label 'regular'

    input:
    file ref from reference

    output:
    file "*.fai" into fai_ch
    file "*.fai" into fai_pychopper_ch

    script:
    """
    samtools faidx $ref
    """
}


process pychopper {

    label "big_mem"

       
    input:
    file fastq from raw_fq


    output:
    file "*" into pychopper_out
    file "pychopper*" into multiQC_0
    file "pychop.fq" into pychopper_ch
    file "pychop.fq" into correction_ch

    shell:
    '''
    cdna_classifier.py -t 16 \
        -r "pychopper_report.pdf" \
        -u "pychopper.unclassified.fq" \
        -w "pychopper.rescued.fq" \
        -S "pychopper.stats" \
        -A "pychopper.scores" \
        "!{fastq}" "pychop.fq"

    '''

}


pychopper_out.flatten().subscribe { it.copyTo("./pychopper/" + name + "/" + name + "_" + it.getName()) }


process mapping_index {

    label 'big_mem'

    input:
    file ref from reference

    output:
    file 'ref.mmi' into mapping
    file 'ref.mmi' into mapping_pychopper
    file 'ref.mmi' into mapping_index_out


    shell:
    '''
    minimap2 -t 16 -d ref.mmi '!{ref}'
    '''
}

mapping_index_out.subscribe { it.copyTo("./mapping/" + name + "/" + it.getName()) }


process mapping_pychopper {

    label 'big_mem'

    input:
    file 'basecall' from correction_ch
    file 'index' from mapping_pychopper

    output:
    file "${name}_mapped_sorted.bam" into mapped_sorted_bam_down_pychopper
    file '*.flagstat' into multiQC_1
    file '*.idxstats' into multiQC_2
    file "*" into mapping_out_pychopper
    file '*sorted.bam' into rseqc_bam
    file '*sorted.bam.bai' into rseqc_bai


    script:
    """
     minimap2 -t 16 -ax splice \
        -uf \
        $index \
        "basecall" > '${name}_mapped.bam' \
         
     samtools sort -@ 16 "${name}_mapped.bam" -o '${name}_mapped_s0rted.bam'
     samtools view -b -F 4 '${name}_mapped_s0rted.bam' > '${name}_mapped_sorted.bam'
     samtools index "${name}_mapped_sorted.bam"
     samtools flagstat "${name}_mapped_sorted.bam" > "${name}_mapped.flagstat"
     samtools idxstats "${name}_mapped_sorted.bam" > "${name}_mapped.idxstats"
     """

}

mapping_out_pychopper.flatten().subscribe { it.copyTo("./mapping_pychopper/" + name + "/" + it.getName()) }


process mapping {

     label 'big_mem'

     input:
     file 'basecall' from raw_fq
     file 'index' from mapping

     output:
     file "mapped_sorted.bam" into mapped_sorted_bam
     file 'mapped_sorted.bam.bai' into mapped_sorted_bai
     file "mapped_sorted.bam" into mapped_raw_bam
     file "mapped_sorted.bam.bai" into mapped_raw_bai


     script:
     """ 
      minimap2 -t 16 -ax splice \
         -uf \
         $index \
         "basecall" > 'mapped.bam' \

      samtools sort -@ 16 "mapped.bam" -o 'mapped_sorted.bam'
      samtools index "mapped_sorted.bam"
      """ 

 }

process quality_control {
    
    label 'big_mem'
    
    input:
    file 'mapped.bam' from mapped_sorted_bam
    file 'mapped.bam.bai' from mapped_sorted_bai
    file basecall from raw_txt


    output:
    file '*' into quality_control_out
    file "*pycoqc.json" into multiQC_3

    shell:
    '''
    fix_pyco_summary.py "!{basecall}" "final_pyco.txt"

    pycoQC -f "final_pyco.txt" \
        -a !{'mapped.bam'} \
        -o "./pycoqc.html" \
        -j "./pycoqc.json" \
        --quiet
    '''

}

quality_control_out.flatten().subscribe { it.copyTo("./pycoqc/" + name + "/" + name + "_" + it.getName()) }


process rseqc {

    label 'big_mem'

    input:
    file bam from rseqc_bam
    file bai from rseqc_bai
    file housekeep from housekeep_file
    file bed from bed_file

    output:
    file "*" into rseqc_out
    file "*" into multiQC_4

    script:

    """
    bam_stat.py -i $bam > "bam_stat"
    geneBody_coverage.py -r $housekeep -i $bam  -o "gene_body_coverage"
    junction_annotation.py -i $bam -o "junction_annotation" -r $bed 2> "junction_annotation_multiqc"
    junction_saturation.py -i $bam -r $bed -o "junction_saturation"
    read_GC.py -i $bam -o "read_gc"
    read_NVC.py -i $bam -o "read_nvc"
    read_quality.py -i $bam -o "read_quality"
    read_distribution.py  -i $bam -r $bed > "read_distribution"
    read_duplication.py -i $bam -o "read_duplication"
    """

}

multiQC_0.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" + it.getName()) }
multiQC_1.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" + it.getName()) }
multiQC_2.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" + it.getName()) }
multiQC_3.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" +  it.getName()) }
multiQC_4.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" + it.getName()) }
rseqc_out.flatten().subscribe { it.copyTo("./rseqc/" + name + "_" + it.getName()) }


process downsampling_pychopper {

label 'regular'

    input:

    each d from downsizes_2
    file "mapped.bam" from mapped_sorted_bam_down_pychopper

    output:
    file "*downsampled*.bam" into downsamples_bam_pychopper
    file "*downsampled*.bam.bai" into downsamples_bai_pychopper
    file "*downsampled*.bam" into downsamples_bam_out_pychopper
    file "*downsampled*.bam.bai" into downsamples_bai_out_pychopper
    file "*downsampled*.bam" into feature_counts_pychopper
    file "*downsampled*.bam" into bambu_pychopper
    file "*downsampled*.bam.bai" into bambu_bai_pychopper
    
    
    shell:
    '''
    samtools view -@ 12 -s '!{d}' -b !{"mapped.bam"} > "!{name}_mapped.downsampled_!{d}.bam"
    samtools index -@ 12 "!{name}_mapped.downsampled_!{d}.bam"
    '''
}

downsamples_bam_out_pychopper.flatten().subscribe { it.copyTo("./downsamples_pychopper/" + name + "/" + it.getName()) }
 
process stringtie {

    label 'regular'

    input:
    file bam from downsamples_bam_pychopper.flatten()
    file bai from downsamples_bai_pychopper.flatten()
    file gtf from annotation

    output:
    file "*.gene_abund*" into final_results_pychopper
    file "*transcripts.gtf" into gff_compare_pychopper
    file "*" into stringtie_out_pychopper

    script:
    """
    stringtie $bam \
        -p 16 \
        -L -o "${bam}.transcripts.gtf" \
        -v \
        -G "${gtf}" \
        -A "${bam}.gene_abund.txt" \
        -C "${bam}.cov_refs.gft" \
        -b "${bam}.ballgown"
    
    mv ${bam}.ballgown/e2t.ctab ${bam}_e2t.ctab
    mv ${bam}.ballgown/e_data.ctab ${bam}_e_data.ctab
    mv ${bam}.ballgown/i2t.ctab ${bam}_i2t.ctab
    mv ${bam}.ballgown/i_data.ctab ${bam}_i_data.ctab
    mv ${bam}.ballgown/t_data.ctab ${bam}_t_data.ctab
    yes | rm -r ${bam}.ballgown
    """

}

stringtie_out_pychopper.flatten().subscribe{ it.copyTo("./stringtie/" + name + "/" + it.getName()) }


process gffcompare_stringtie {

    label 'regular'

    input:
    file transcripts from gff_compare_pychopper
    file gtf from annotation

    output:
    file '*' into gff_out_pychopper

    shell:
    '''
    gffcompare -r "!{gtf}" "!{transcripts}"
    mv gffcmp.annotated.gtf !{transcripts}.annotated.gtf
    mv gffcmp.stats !{transcripts}.stats
    mv gffcmp.tracking !{transcripts}.tracking
    mv gffcmp.loci !{transcripts}.loci
    mv "gffcmp.!{transcripts}.tmap" "!{transcripts}.tmap"
    mv "gffcmp.!{transcripts}.refmap" "!{transcripts}.refmap"
    '''

}

gff_out_pychopper.flatten().subscribe { it.copyTo("./gffcompare_stringtie/" + name + "/" + it.getName()) }


process feature_counts {

    label 'regular'

    input:
    file bam from feature_counts_pychopper
    file gtf from annotation

    output:
    file "*" into feature_counts_out_pychopper
    file "*counts.txt" into feature_counts_final_pychopper

    script:
    """
    featureCounts \
     -L \
     -t exon \
     -g gene_id \
     -a $gtf \
     --extraAttributes 'gene_name' \
     -o "${bam}_counts.txt" $bam
    """

}


feature_counts_out_pychopper.flatten().subscribe{ it.copyTo("./feature_counts/" + name + "/" + it.getName()) }



process bambu {
 
    label 'big_mem'
 
    input:
    file gtf from annotation
    file ref from reference
    file bam from mapped_raw_bam
    file ref_fai from fai_pychopper_ch
    file bai from mapped_raw_bai

    output:
    file "no_discovery/genes/*" into bambu_out_pychopper_1
    file "no_discovery/transcripts/*" into bambu_out_pychopper_2
    file "no_discovery/annotations/*" into bambu_out_pychopper_3
    file "novel/genes/*" into bambu_out_pychopper_4
    file "novel/transcripts/*" into bambu_out_pychopper_5
    file "novel/annotations/*" into bambu_out_pychopper_6
    file "no_discovery/genes/*" into bambu_final_results_pychopper
    file "novel/annotations/*" into bambu_gff_pychopper_ch 
 
    shell:
    '''
    mkdir -p no_discovery
    mkdir -p novel
 
    bambu_script.R "!{bam}" "!{ref}" "!{gtf}" "!{name}"
    bambu_script_novel.R "!{bam}" "!{ref}" "!{gtf}" "!{name}"

    mkdir -p no_discovery/genes/
    mkdir -p no_discovery/transcripts/
    mkdir -p no_discovery/annotations/
    mkdir -p novel/genes/
    mkdir -p novel/transcripts/
    mkdir -p novel/annotations/

    mv no_discovery/counts_g* no_discovery/genes/"!{bam}_counts_gene.txt"
    mv no_discovery/counts_t* no_discovery/transcripts/"!{bam}_counts_transcript.txt"
    mv no_discovery/extended* no_discovery/annotations/"!{bam}_extended_annotations.gtf"
 
 
    mv novel/counts_g* novel/genes/"!{bam}_counts_gene.txt"
    mv novel/counts_t* novel/transcripts/"!{bam}_counts_transcript.txt"
    mv novel/extended* novel/annotations/"!{bam}_extended_annotations.gtf"
 
   '''
}
 
bambu_out_pychopper_1.flatten().subscribe{ it.copyTo("./bambu/" + name + "/no_discovery/genes/" + it.getName()) }
bambu_out_pychopper_2.flatten().subscribe{ it.copyTo("./bambu/" + name + "/no_discovery/transcripts/" + it.getName()) }
bambu_out_pychopper_3.flatten().subscribe{ it.copyTo("./bambu/" + name + "/no_discovery/annotations/" + it.getName()) }
bambu_out_pychopper_4.flatten().subscribe{ it.copyTo("./bambu/" + name + "/novel/genes/" + it.getName()) }
bambu_out_pychopper_5.flatten().subscribe{ it.copyTo("./bambu/" + name + "/novel/transcripts/" + it.getName()) }
bambu_out_pychopper_6.flatten().subscribe{ it.copyTo("./bambu/" + name + "/novel/annotations/" + it.getName()) }



process gffcompare_bambu {
 
    label 'regular'

    input:
    file transcripts from bambu_gff_pychopper_ch.flatten()
    file gtf from annotation
 
    output:
    file '*' into gff_bambu_pychopper_out
 
    shell:
    '''
    gffcompare -r "!{gtf}" "!{transcripts}"
    mv gffcmp.annotated.gtf !{transcripts}.annotated.gtf
    mv gffcmp.stats !{transcripts}.stats
    mv gffcmp.tracking !{transcripts}.tracking
    mv gffcmp.loci !{transcripts}.loci
    mv "gffcmp.!{transcripts}.tmap" "!{transcripts}.tmap"
    mv "gffcmp.!{transcripts}.refmap" "!{transcripts}.refmap"
    '''
 
}
 
gff_bambu_pychopper_out.flatten().subscribe { it.copyTo("./gffcompare_bambu/" + name + "/" + it.getName()) }


process collector {

    
    label 'regular'
    
    input:
    file four from final_results_pychopper.toSortedList({ a,b -> a.baseName <=> b.baseName  })
    file five from bambu_final_results_pychopper.toSortedList({ a,b -> a.baseName <=> b.baseName  })
    file six from feature_counts_final_pychopper.toSortedList({ a,b -> a.baseName <=> b.baseName  })
    
    output:
    file four into final_results_pychopper_col
    file five into bambu_final_results_pychopper_col
    file six into feature_counts_final_pychopper_col
    

    shell:
    """
    echo "hello"
    """
}


process final_results {

    label 'regular'

    input:
    file stringtie from final_results_pychopper_col.flatten()
    file bambu from bambu_final_results_pychopper_col.flatten()
    file features from feature_counts_final_pychopper_col.flatten()
    
    output:
    file '*final_results' into final_out_pychopper


    shell:
    '''
    final_results_long.py "!{bambu}" "!{stringtie}" "!{features}" "!{stringtie}.final_results"
    '''

}

final_out_pychopper.flatten().subscribe{ it.copyTo("./final_results/" + name + "/" + it.getName()) }
