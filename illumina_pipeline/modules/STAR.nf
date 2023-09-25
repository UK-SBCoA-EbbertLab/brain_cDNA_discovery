process MAKE_STAR_INDEX {

    label 'huge'

    input:
        path(reference)
        path(annotation)
        val(overhang)

    output:
        path "./STAR_indexes"


    script:
        """
        STAR --runThreadN 32 \
            --runMode genomeGenerate \
            --genomeDir ./STAR_indexes \
            --genomeFastaFiles $reference \
            --sjdbGTFfile $annotation \
            --sjdbOverhang $overhang \
 
        """
}

process STAR_MAPPING {

    publishDir "results/${params.out_dir}/STAR/", mode: 'copy', overwrite: true

    label 'medium_large'

    input:
        val(id)
        path(trimmed_R1)
        path(trimmed_R2)
        path(index)

    output:
        val "$id", emit: id
        path "${id}_Aligned.toTranscriptome.out.bam", emit: bam_trans
        path "*", emit: outty

    script:
        """
        STAR --genomeDir "${index}" \
            --runThreadN 16 \
            --readFilesIn "${trimmed_R1}" "${trimmed_R2}" \
            --quantMode TranscriptomeSAM \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix ./"${id}_"

        samtools flagstat "${id}_Aligned.toTranscriptome.out.bam" > "${id}_Aligned.toTranscriptome.out.bam.flagstat"
        """
}
