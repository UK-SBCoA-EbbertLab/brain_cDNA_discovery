#!/usr/bin/env nextflow

params.file_1 = 'raw_data/Brain-421Wet_S2_R1_001.fastq.gz'
params.file_2 = 'raw_data/Brain-421Wet_S2_R2_001.fastq.gz'
params.ref = 'references/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.gtf = 'references/Homo_sapiens.GRCh38.104.gtf'
params.name = "sample"
params.bed = "references/hg38_ncbi_08_29_2021.bed"
params.housekeep = "references/hg38.HouseKeepingGenes.bed"

file_1 = file(params.file_1)
file_2 = file(params.file_2)
ref = file(params.ref)
gtf = file(params.gtf)
name = params.name



downsizes = Channel.from(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
annotation = Channel.value(gtf)
reference = Channel.value(ref)
file_1_ch = Channel.fromPath(file_1)
file_2_ch = Channel.fromPath(file_2)
bed_file = Channel.from(file(params.bed))
housekeep_file = Channel.from(file(params.housekeep))

process trim_galore {
 
    label "regular"

    input:
    file file_R1 from file_1_ch
    file file_R2 from file_2_ch

    output:
    file "*" into trimmed_out
    file "*val_1*.fq" into trimmed_ch_1
    file "*val_2*.fq" into trimmed_ch_2
    file "*report.txt" into multiQC_0
    file "*fastqc.zip" into multiQC_1

    shell:
    '''
    mv "!{file_R1}" "!{name}_S2_R1_001.fastq.gz"
    mv "!{file_R2}" "!{name}_S2_R2_001.fastq.gz"

    trim_galore --paired --illumina -j 8 --fastqc "!{name}_S2_R1_001.fastq.gz" "!{name}_S2_R2_001.fastq.gz" -o ./ --basename "!{name}"
    pigz -p16 -d "!{name}_val_1.fq.gz"
    pigz -p16 -d "!{name}_val_2.fq.gz"
    '''

}

trimmed_out.flatten().subscribe { it.copyTo("./trim_galore/" + name + "/" + it.getName()) }


process make_index {

    label "big_mem"

    input:
    file ref from reference
    file gtf from annotation

    output:
    path "./STAR_indexes" into index_ch
    path "./STAR_indexes" into index_out

    script:
    """
    STAR --runThreadN 32 \
         --runMode genomeGenerate \
         --genomeDir ./STAR_indexes \
         --genomeFastaFiles $ref \
         --sjdbGTFfile $gtf \
         --sjdbOverhang 149 \          
    """

}

process mapping {

    label "big_mem"

    input:
    file file_1 from trimmed_ch_1
    file file_2 from trimmed_ch_2
    path index from index_ch

    output:
    file "*" into mapping_out
    file '*sorted.bam' into mapping_ch
    file "*Log.final.out" into multiQC_2
    file "*sorted.bam.bai" into rseqc_bai
    file '*sorted.bam' into rseqc_bam
    file "*flagstat" into multiQC_3
    file "*idxstats" into multiQC_4

    script:
    """
    STAR --genomeDir $index \
         --runThreadN 16 \
         --readFilesIn $file_1 $file_2 \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ./"${name}_"

    samtools sort -@ 16 "${name}_Aligned.sortedByCoord.out.bam" -o '${name}_mapped_sorted.bam'
    samtools index '${name}_mapped_sorted.bam'
    samtools flagstat '${name}_mapped_sorted.bam' > "${name}_mapped_sorted.flagstat"
    samtools idxstats "${name}_mapped_sorted.bam" > "${name}_mapped_sorted.idxstats"

    yes | rm -r "${name}_Aligned.sortedByCoord.out.bam"
    """

}

mapping_out.flatten().subscribe { it.copyTo("./star/" + name + "/" + it.getName()) }

process rseqc {

    label 'regular'

    input:
    file bam from rseqc_bam
    file housekeep from housekeep_file
    file bed from bed_file
    file bai from rseqc_bai

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
    inner_distance.py -i $bam -o "inner_distance" -r $bed
    read_NVC.py -i $bam -o "read_nvc"
    read_quality.py -i $bam -o "read_quality"
    read_distribution.py  -i $bam -r $bed > "read_distribution"
    read_duplication.py -i $bam -o "read_duplication"
    """

}

multiQC_0.flatten().subscribe { it.copyTo("./multiqc/" + it.getName()) }
multiQC_1.flatten().subscribe { it.copyTo("./multiqc/" + it.getName()) }
multiQC_2.flatten().subscribe { it.copyTo("./multiqc/" + it.getName()) }
multiQC_3.flatten().subscribe { it.copyTo("./multiqc/" + it.getName()) }
multiQC_4.flatten().subscribe { it.copyTo("./multiqc/" + it.getName()) }
multiQC_5.flatten().subscribe { it.copyTo("./multiqc/" + name + "_" + it.getName()) }
rseqc_out.flatten().subscribe { it.copyTo("./rseqc/" + name + "_" + it.getName()) }

process downsampling {

    label "regular"

    input:
    each d from downsizes
    file mapped from mapping_ch

    output:
    file "*downsampled*.bam" into downsampled_bam_ch
    file "*downsampled*.bam.bai" into downsampled_bai_ch
    file "*downsampled*.bam" into feature_counts
    file "*" into downsampled_out


    script:
    """
    samtools view -@ 12 -s '${d}' -b $mapped > "${name}_mapped.downsampled_${d}.bam"
    samtools index -@ 12 "${name}_mapped.downsampled_${d}.bam"
    """

}

downsampled_out.flatten().subscribe { it.copyTo("./downsamples/" + name + "/" + it.getName()) }


process stringtie {

    label "regular"

    input:
    file bam from downsampled_bam_ch.flatten()
    file bai from downsampled_bai_ch.flatten()
    file gtf from annotation

    output:
    file "*.gene_abund*" into final_results
    file "*transcripts.gtf" into gffcompare
    file "*" into stringtie_out

    shell:
    '''
    stringtie "!{bam}" \
        -p 16 \
        -o "!{bam}_transcripts.gtf" \
        -v \
        -G "!{gtf}" \
        -A "!{bam}.gene_abund.txt" \
        -C "!{bam}.cov_refs.gtf" \
        -b "!{bam}.ballgown"
    
    mv !{bam}.ballgown/e2t.ctab !{bam}_e2t.ctab
    mv !{bam}.ballgown/e_data.ctab !{bam}_e_data.ctab
    mv !{bam}.ballgown/i2t.ctab !{bam}_i2t.ctab
    mv !{bam}.ballgown/i_data.ctab !{bam}_i_data.ctab
    mv !{bam}.ballgown/t_data.ctab !{bam}_t_data.ctab
    yes | rm -r !{bam}.ballgown
    '''
}

stringtie_out.flatten().subscribe{ it.copyTo("./stringtie/" + name + "/" + it.getName()) }


process gffcompare {

    label "regular"

    input:
    file transcripts from gffcompare.flatten()
    file gtf from annotation

    output:
    file "*" into gff_out

    shell:
    '''
    gffcompare -r "!{gtf}" "!{transcripts}"
    mv gffcmp.annotated.gtf !{transcripts}.annotated.gtf
    mv gffcmp.stats !{transcripts}.stats
    mv gffcmp.tracking !{transcripts}.tracking
    mv gffcmp.loci !{transcripts}.loci
    '''
}

gff_out.flatten().subscribe { it.copyTo("./gffcompare_stringtie/" + name + "/" + it.getName()) }


process feature_counts {


    label "regular"

    input:
    file mapped from feature_counts
    file gtf from annotation

    output:
    file "*feature_counts.txt" into feature_counts_ch
    file "*" into feature_counts_out

    script:
    """
    featureCounts \
    -T 12 \
    -p \
    -t exon \
    -g gene_id \
    -a $gtf \
    --extraAttributes 'gene_name' \
    -o "${mapped}_feature_counts.txt" "${mapped}"
    """
}

feature_counts_out.flatten().subscribe { it.copyTo("./feature_counts/" + name + "/" + it.getName()) }


process collector {

    label 'regular'

    input:
    file one from final_results.toSortedList({ a,b -> a.baseName <=> b.baseName  })
    file three from feature_counts_ch.toSortedList({ a,b -> a.baseName <=> b.baseName  })

    output:
    file one into final_results_col
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
    file features from feature_counts_final_col.flatten()
    

    output:
    file '*final_results' into final_out
    
    shell:
    '''
    final_results_short.py "!{stringtie}" "!{features}" "!{stringtie}.final_results"
    '''
}

final_out.flatten().subscribe{ it.copyTo("./final_results/" + name + "/" + it.getName()) }

