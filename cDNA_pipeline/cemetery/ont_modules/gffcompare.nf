process GFFCOMPARE_STRINGTIE {

    publishDir 'results/gffcompare_stringtie/', mode: 'copy', overwrite: false

    label 'small'

    input:
        val(id)
        path(transcripts)
        path(gtf)

    output:
        path "*"

    script:
    """
    gffcompare -r $gtf $transcripts -o "${id}_stringtie"
    """
}


process GFFCOMPARE_BAMBU {

    publishDir 'results/gffcompare_bambu/', mode: 'copy', overwrite: false

    label 'small'

    input:
        val(id)
        path(transcripts)
        path(gtf)

    output:
        path "*"

    script:
    """
    gffcompare -r $gtf $transcripts -o "${id}_bambu"
    """
}

