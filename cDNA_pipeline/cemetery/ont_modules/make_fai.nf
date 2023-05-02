process MAKE_FAI {

    publishDir 'results/make_fai/', mode: 'copy', overwrite: false

    label 'tiny'

    input:
        path ref

    output:
        path '*.fai'

    script:
    """
    samtools faidx $ref
    """
}




