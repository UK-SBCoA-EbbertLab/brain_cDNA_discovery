process MAKE_TRANSCRIPTOME {

    publishDir "results/${params.out_dir}/transcriptome/", mode: "copy", overwrite: true

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
