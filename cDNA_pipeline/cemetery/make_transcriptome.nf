process MAKE_TRANSCRIPTOME {

    publishDir "results/${params.out_dir}/transcriptome/"

    label 'medium'

    input:
        path genome
        path genome_index
        path annotation

    output:
        path('transcriptome.fa')

    script:
        """
        gffread -w transcriptome.fa -g $genome $annotation
        """
}
