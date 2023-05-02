process MAKE_INDEX_cDNA {

    label 'small'

    input:
        path ref

    output:
        path 'ref.mmi'


    script:
        """
        minimap2 -t 8 -d ref.mmi $ref
        """
}
