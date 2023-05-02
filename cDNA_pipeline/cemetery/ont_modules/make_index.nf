process MAKE_INDEX {

    publishDir 'results/mapping/', mode: 'copy', overwrite: false

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
