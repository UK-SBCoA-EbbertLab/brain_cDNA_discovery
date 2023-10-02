process FILTER_UNIQUE_BAM {

    publishDir "results/${params.out_dir}/unique_alignment_bams/", mode: 'copy', overwrite: true

    label "medium"

    input:
        tuple val(id), path(bam)

    output:
        val "$id", emit: id
        path "${id}_transcriptome_uniquely_aligned_mapq_255.bam", emit: bam

    script:
        """
        
        samtools view -b -h -q 255 "${bam}" > "${id}_transcriptome_uniquely_aligned_mapq_255.bam"

        """

}
