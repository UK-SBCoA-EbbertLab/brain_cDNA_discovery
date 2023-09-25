process STRINGTIE_ONT_cDNA {

    publishDir "results/${params.out_dir}/stringtie_nanopore_cDNA/" , mode: 'copy', overwrite: true

    label 'medium'

    input:
        val(id)
        path(bam)
        path(bai)
        path(gtf)

    output:
        path "*"

    script:
        """

        mkdir "${id}"

        stringtie $bam \
            -p 16 \
            -L -o "./${id}/${id}.transcripts.gtf" \
            -v \
            -e \
            -G "${gtf}" \
            -A "./${id}/${id}.gene_abund.txt" \
            -C "./${id}/${id}.cov_refs.gtf" 
    """

}
