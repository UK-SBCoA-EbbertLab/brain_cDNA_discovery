process FILTER_UNIQUE_BAM {

    publishDir "results/${params.out_dir}/unique_alignment_bams/", mode: 'copy', overwrite: true

    label "medium_large"

    input:
        tuple val(id), path(bam_file)

    output:
        val "$id", emit: id
        path "${id}_transcriptome_uniquely_aligned_mapq_255.bam", emit: bam
        path "*", emit: outty

    script:
        """
        
        samtools view -b -h -q 255 "${bam_file}" > "${id}_transcriptome_uniquely_aligned_mapq_255.bam"

        samtools flagstat "${id}_transcriptome_uniquely_aligned_mapq_255.bam" > "${id}_transcriptome_uniquely_aligned_mapq_255.bam.flagstat"
        """

}
