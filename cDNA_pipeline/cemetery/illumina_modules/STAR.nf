process MAKE_STAR_INDEX {

    publishDir 'results/mapping/', mode: 'copy', overwrite: false

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

    publishDir 'results/mapping/', mode: 'copy', overwrite: false

    label 'medium_large'

    input:
        val(id)
        path(trimmed_R1)
        path(trimmed_R2)
        path(index)

    output:
        val "$id", emit: id
        path "${id}_STAR_aligned_transcriptome_mapped.bam", emit: bam
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

    mv '${id}_Aligned.toTranscriptome.out.bam' '${id}_STAR_aligned_transcriptome_all.bam'
    samtools view -b -F 4 -@ 12 '${id}_STAR_aligned_transcriptome_all.bam' > '${id}_STAR_aligned_transcriptome_mapped.bam'
    samtools flagstat "${id}_STAR_aligned_transcriptome_mapped.bam" > "${id}_STAR_aligned_transcriptome_mapped.flagstat"
    """
}
