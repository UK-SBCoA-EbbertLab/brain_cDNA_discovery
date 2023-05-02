process STRINGTIE_DISCOVERY {

    publishDir 'results/stringtie/', mode: 'copy', overwrite: false

    label 'medium_small'

    input:
        val(id)
        path(bam)
        path(index)
        path(gtf)

    output:
        val(id), emit: id
        path "*transcripts.gtf", emit: annotation
        path "*", emit: out

    script:
    """
    stringtie $bam \
        -p 12 \
        -L -o "${id}.transcripts.gtf" \
        -v \
        -G "${gtf}" \
        -A "${id}.gene_abund.txt" \
        -C "${id}.cov_refs.gft" \
        -b "${id}.ballgown"
    """

}

