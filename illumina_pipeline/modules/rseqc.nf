process RSEQC {

    publishDir "results/${params.out_dir}/QC/RseQC/", mode: 'copy', overwrite: true
    
    label "large"

    input:
        path(bam_illumina)
        path(bai_illumina)
        path(bam_nanopore)
        path(bai_nanopore)
        path(housekeeping)

    output:
        path "*", emit: multiQC

    shell:
        '''
        concatenated_bams="!{bam_illumina} !{bam_nanopore}"

        final_bams="$(tr ' ', ',' <<<$concatenated_bams)"

        geneBody_coverage.py -i $final_bams -r "!{housekeeping}" -o "geneBody_coverage.py"
        '''
}

process RSEQC_ONT_ONLY {

    publishDir "results/${params.out_dir}/QC/RseQC/", mode: "copy", overwrite: true
 
    label "large"

    input:
        path(bam_nanopore)
        path(bai_nanopore)
        path(housekeeping)

    output:
        path "*", emit: multiQC

    shell:
        '''
        concatenated_bams="!{bam_nanopore}"

        final_bams="$(tr ' ', ',' <<<$concatenated_bams)"

        geneBody_coverage.py -i $final_bams -r "!{housekeeping}" -o "geneBody_coverage.py"
        '''
}
