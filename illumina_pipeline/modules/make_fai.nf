process MAKE_FAI {

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
