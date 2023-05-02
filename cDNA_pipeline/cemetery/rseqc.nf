process RSEQC {

    publishDir "results/${params.out_dir}/RseQC/"
 
    label "large"

    input:
        path(bam_illumina)
        path(bam_nanopore)
        path(bai_nanopore)
        path(housekeeping)

    output:
        path "*"

    shell:
        '''
        samtools index !bam_illumina

        concatenated_bams="!{bam_illumina} !{bam_nanopore}"

        final_bams="$(tr ' ', ',' <<<$concatenated_bams)"

        geneBody_coverage.py -i $final_bams -r !housekeeping -o "geneBody_coverage.py"

        '''
}
