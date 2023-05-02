process MAKE_STAR_INDEX {

    label 'huge'

    input:
        path reference
        path annotation

    output:
        path "./STAR_indexes"


    script:
        """
        STAR --runThreadN 32 \
             --runMode genomeGenerate \
            --genomeDir ./STAR_indexes \
            --genomeFastaFiles $reference \
            --sjdbGTFfile $annotation \
            --sjdbOverhang 149 \
 
        """
}

process STAR_MAPPING_TRANSCRIPTOME {

    publishDir "results/${params.out_dir}/mapping/"

    label 'medium_large'

    input:
        val(id)
        path(trimmed_R1)
        path(trimmed_R2)
        path(index)

    output:
        val "$id", emit: id
        path "${id}_mapped_filtered_trans.bam", emit: bam_trans
        path "${id}_mapped_filtered_gen.bam", emit: bam_gen
        path "*Log.final.out", emit: QC_3
        path "*stat", emit: QC_4

    script:
        """
        STAR --genomeDir "${index}" \
            --runThreadN 16 \
            --quantMode TranscriptomeSAM \
            --readFilesIn "${trimmed_R1}" "${trimmed_R2}" \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix ./"${id}_"


        mv "${id}_Aligned.out.bam" '${id}_all_gen.bam'
        samtools flagstat '${id}_all_gen.bam' > '${id}_all_gen.flagstat'
        samtools view -b -q 10 -F 2304 -@ 12  '${id}_all_gen.bam' > '${id}_mapped_filtered_gen.bam'
        samtools flagstat "${id}_mapped_filtered_gen.bam" > "${id}_mapped_filtered_gen.flagstat"


        mv '${id}_Aligned.toTranscriptome.out.bam' '${id}_all_trans.bam'
        samtools flagstat '${id}_all_trans.bam' > '${id}_all_trans.flagstat'
        samtools view -b -q 10 -F 2304 -@ 12  '${id}_all_trans.bam' > '${id}_mapped_filtered_trans.bam'
        samtools flagstat "${id}_mapped_filtered_trans.bam" > "${id}_mapped_filtered_trans.flagstat"
        """
}
