#!/usr/bin/env nextflow

params.bam = 'raw_data/m64014_190506_005857.subreads.bam'
params.pbi = 'raw_data/m64014_190506_005857.subreads.bam.pbi'
params.adapters = 'raw_data/m64014_190506_005857.adapters.fasta'

params.ref = 'references/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.gtf = 'references/Homo_sapiens.GRCh38.104.gtf'
params.name = "sample"
params.bed = "references/hg38_ncbi_08_29_2021.bed"
params.housekeep = "references/hg38.HouseKeepingGenes.bed"


 
downsizes = Channel.from(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
annotation = Channel.value(file(params.gtf))
reference = Channel.value(file(params.ref))
adapters = Channel.value(file(params.adapters))
rawdata = Channel.value(file(params.bam))
pbi = Channel.value(file(params.pbi))
name = params.name
bed_file = Channel.from(file(params.bed))
housekeep_file = Channel.from(file(params.housekeep))


process ccs {

    label "regular"


    input:
    each x from 1..50
    file bam from rawdata
    file index from pbi

    output:
    file "*bam" into gather_ccs_ch
    file "*report*" into gather_ccs_report_ch
    file "*pbi" into gather_pbi_ch

    script:
    """
    ccs $bam "${name}_chunk_${x}.bam" --chunk "${x}"/50  --min-passes 1 --min-rq 0.9  -j 12 --report-file "${name}_ccs_report_chunk_${x}"
    """


}


process gather_chunks {

    label "regular"

    input:
    file "chunk*" from gather_ccs_ch.collect()
    file "pbi_chunk*" from gather_pbi_ch.collect()
    file "report*" from gather_ccs_report_ch.collect()

    output:
    file "${name}.ccs.bam" into lima_bam_ch
    file "${name}.ccs.bam.pbi" into lima_pbi_ch
    file "*" into gather_chunks_out


    script:
    """
    samtools merge -@ 12 "${name}.ccs.bam" chunk*
    pbindex "${name}.ccs.bam"
    combine_ccs_report.py report* > "${name}.ccs_report.txt"
    """
}

gather_chunks_out.flatten().subscribe { it.copyTo("./ccs/" + name + "/" + it.getName()) }


process make_fai {

    label 'regular'
 
    input:
    file ref from reference
 
    output:
    file "*.fai" into fai_ch
 
    script:
    """
    samtools faidx $ref
    """
}


process lima {

    label "regular"

    input:
    file bam from lima_bam_ch
    file pbi from lima_pbi_ch
    file primers from adapters

    output:
    file "*bam" into refine_bam_ch
    file "*pbi" into refine_pbi_ch
    file "*" into lima_out_ch
    file "*summary" into multiQC_1
    file "*counts" into multiQC_2

    script:
    """
    lima --num-threads 12 --dump-clips --peek-guess --isoseq $bam $primers "demux.bam"
    """

}

lima_out_ch.flatten().subscribe { it.copyTo("./lima/" + name + "/" + name + "_" + it.getName()) }

process refine {

    label "regular"

    input:
    file bam from refine_bam_ch
    file pbi from refine_pbi_ch
    file primers from adapters

    output:
    file "*bam" into cluster_bam_ch
    file "*pbi" into cluster_pbi_ch
    file "*fastq" into mapping_fastq_ch
    file "*" into refine_out



    script:
    """
    isoseq3 refine $bam $primers "flnc.bam" --require-polya
    samtools fastq -@ 12 "flnc.bam" > "flnc.fastq"

    """

}

refine_out.flatten().subscribe { it.copyTo("./refine/" + name + "/" + name + "_" + it.getName()) }

process cluster {

    label "regular"

    input:
    file bam from cluster_bam_ch
    file pbi from cluster_pbi_ch
    
    output:
    file "*" into cluster_out

    script:
    """
    
    isoseq3 cluster $bam "unpolished.bam" --use-qvs
    gzip -d *.gz
    
    """

}

cluster_out.flatten().subscribe { it.copyTo("./cluster/" + name + "/" + name + "_" + it.getName()) }


process mapping_index {
 
    label 'big_mem'

    input:
    file ref from reference
    
    output:
    file 'ref.mmi' into mapping_ch
    file 'ref.mmi' into mapping_index_out

    shell:
    '''
    minimap2 -t 16 -d ref.mmi '!{ref}'
    '''
}

mapping_index_out.subscribe { it.copyTo("./mapping/" + name + "/" + it.getName()) }

process mapping {

    label "big_mem"

    input:
    file 'fastq' from mapping_fastq_ch
    file 'index' from mapping_ch

    output:
    file "*mapped_sorted.bam" into mapped_sorted_bam
    file "*mapped_sorted.bam" into mapped_sorted_bam_down
    file '*mapped_sorted.bam.bai' into mapped_sorted_bai
    file "*mapped_sorted.bam" into quality_control
    file "*" into mapping_out 
    file "*mapped_sorted.bam" into rseqc_bam
    file "*mapped_sorted.bam.bai" into rseqc_bai
    file "*flagstat" into multiQC_3
    file "*idxstats" into multiQC_4

    script:
    """
    echo "${name}"
    minimap2 -t 16 -ax splice \
        -uf \
        $index \
        $fastq > "${name}_mapped.bam" \

    samtools sort -@ 12 "${name}_mapped.bam"  -o '${name}_mapped_sorted.bam'
    samtools index "${name}_mapped_sorted.bam"
    samtools flagstat "${name}_mapped_sorted.bam" > "${name}_mapped.flagstat"
    samtools idxstats "${name}_mapped_sorted.bam" > "${name}_mapped.idxstats"
    """
}

mapping_out.flatten().subscribe { it.copyTo("./mapping/" + name + "/" + name + "_" + it.getName()) }

process rseqc {

    label 'regular'

    input:
    file bam from rseqc_bam
    file bai from rseqc_bai
    file housekeep from housekeep_file
    file bed from bed_file

    output:
    file "*" into rseqc_out
    file "*" into multiQC_5

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

multiQC_1.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" + it.getName()) }
multiQC_2.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" + it.getName()) }
multiQC_3.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" +  it.getName()) }
multiQC_4.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" + it.getName()) }
multiQC_5.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" + it.getName()) }
rseqc_out.flatten().subscribe { it.copyTo("./rseqc/" + name + "_" + it.getName()) }



process downsampling {

    label 'regular'

    input:
    each d from downsizes
    file bam from mapped_sorted_bam_down

    output:
    file "*downsampled*.bam" into downsamples_bam
    file "*downsampled*.bam.bai" into downsamples_bai
    file "*downsampled*.bam" into feature_counts_ch
    file "*" into downsampling_out
    file "*downsampled*.bam" into bambu
    file "*downsampled*.bam.bai" into bambu_bai
  
    script:
    """
    samtools view -@ 12 -s '${d}' -b $bam > "mapped.downsampled_${d}.bam"
    samtools index -@ 12 "mapped.downsampled_${d}.bam"
    """
}

downsampling_out.flatten().subscribe { it.copyTo("./downsampling/" + name + "/" + name + "_" + it.getName()) }

process stringtie {

label 'regular'

    input:
    file bam from downsamples_bam.flatten()
    file bai from downsamples_bai.flatten()
    file gtf from annotation

    output:
    file "*.gene_abund*" into final_results
    file "*transcripts.gtf" into gff_compare
    file "*" into stringtie_out

    shell:
    '''
    stringtie "!{bam}" \
    -p 12 \
    -L -o "!{bam}.transcripts.gtf" \
    -v \
    -G "!{gtf}" \
    -A "!{bam}.gene_abund.txt" \
    -C "!{bam}.cov_refs.gft" \
    -b "!{bam}.ballgown"

    mv "!{bam}.ballgown/e2t.ctab" "!{bam}_e2t.ctab"
    mv "!{bam}.ballgown/e_data.ctab" "!{bam}_e_data.ctab"
    mv "!{bam}.ballgown/i2t.ctab" "!{bam}_i2t.ctab"
    mv "!{bam}.ballgown/i_data.ctab" "!{bam}_i_data.ctab"
    mv "!{bam}.ballgown/t_data.ctab" "!{bam}_t_data.ctab"
    yes | rm -r "!{bam}.ballgown"
    '''
}

stringtie_out.flatten().subscribe { it.copyTo("./stringtie/" + name + "/" + name + "_" + it.getName()) }



process gffcompare_stringtie {

    label 'regular'

    input:
    file transcripts from gff_compare.flatten()
    file gtf from annotation

    output:
    file '*' into gff_out

    shell:
    '''
    gffcompare -r "!{gtf}" "!{transcripts}"
    mv gffcmp.annotated.gtf !{transcripts}.annotated.gtf
    mv gffcmp.stats !{transcripts}.stats
    mv gffcmp.tracking !{transcripts}.tracking
    mv gffcmp.loci !{transcripts}.loci
    '''
}

gff_out.flatten().subscribe { it.copyTo("./gffcompare_stringtie/" + name + "/" + name + "_" + it.getName()) }

/*
process feature_counts {

    label 'regular'

    input:
    file bam from feature_counts_ch
    file gtf from annotation

    output:
    file "*_counts.txt" into final_feature_ch
    file "*" into feature_counts_out

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

feature_counts_out.flatten().subscribe { it.copyTo("./feature_counts/" + name + "/" + name + "_" + it.getName()) }


process bambu {
 
    label 'big_mem'
 
    input:
    file gtf from annotation
    file ref from reference
    file bam from bambu
    file ref_fai from fai_ch
    file bai from bambu_bai

    output:
    file "no_discovery/genes/*" into bambu_out_1
    file "no_discovery/transcripts/*" into bambu_out_2
    file "no_discovery/annotations/*" into bambu_out_3
    file "novel/genes/*" into bambu_out_4
    file "novel/transcripts/*" into bambu_out_5
    file "novel/annotations/*" into bambu_out_6
    file "no_discovery/genes/*" into bambu_final_results
    file "novel/annotations/*" into bambu_gff_ch 
 
    shell:
    '''
    mkdir -p novel
    mkdir -p no_discovery
    
    bambu_script.R "!{bam}" "!{ref}" "!{gtf}" "!{name}"
 
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

bambu_out_1.flatten().subscribe{ it.copyTo("./bambu/" + name + "/no_discovery/genes/" + name + "_" + it.getName()) }
bambu_out_2.flatten().subscribe{ it.copyTo("./bambu/" + name + "/no_discovery/transcripts/" + name + "_" + it.getName()) }
bambu_out_3.flatten().subscribe{ it.copyTo("./bambu/" + name + "/no_discovery/annotations/" + name + "_" + it.getName()) }
bambu_out_4.flatten().subscribe{ it.copyTo("./bambu/" + name + "/novel/genes/" + name + "_" + it.getName()) }
bambu_out_5.flatten().subscribe{ it.copyTo("./bambu/" + name + "/novel/transcripts/" + name + "_" + it.getName()) }
bambu_out_6.flatten().subscribe{ it.copyTo("./bambu/" + name + "/novel/annotations/" + name + "_" + it.getName()) }


process gffcompare_bambu {
 
    label 'regular'

    input:
    file transcripts from bambu_gff_ch.flatten()
    file gtf from annotation
 
    output:
    file '*' into gff_bambu_out
 
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
 
gff_bambu_out.flatten().subscribe { it.copyTo("./gffcompare_bambu/" + name + "/" + name + "_" + it.getName()) }


process collector {

    label 'regular'

    input:
    file one from final_results.toSortedList({ a,b -> a.baseName <=> b.baseName  })
    file two from bambu_final_results.toSortedList({ a,b -> a.baseName <=> b.baseName  })
    file three from final_feature_ch.toSortedList({ a,b -> a.baseName <=> b.baseName  })

    output:
    file one into final_results_col
    file two into bambu_final_results_col
    file three into feature_counts_final_col


    shell:
    """
    echo "hello"
    """
}

process final_results {

    label 'regular'

    input:
    file stringtie from final_results_col.flatten()
    file bambu from bambu_final_results_col.flatten()
    file features from feature_counts_final_col.flatten()


    output:
    file '*final_results' into final_out

    shell:
    '''
    final_results_long.py "!{bambu}" "!{stringtie}" "!{features}" "!{stringtie}.final_results"

    echo "final_results_long.py !{bambu} !{stringtie} !{features} !{stringtie}.final_results" > test_final.txt
    '''
}

final_out.flatten().subscribe{ it.copyTo("./final_results/" + name + "/" + name + "_" + it.getName()) }
*/
