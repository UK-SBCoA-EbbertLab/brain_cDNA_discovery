process MULTIQC {

    publishDir "results/${params.out_dir}/multiqc/"

    label 'tiny'

    input:
        path(QC_1)
        path(QC_2)
        path(QC_3)
        path(QC_4)
        path(multiqc_config)
    
    output: 
        path "*"

    script:
        """    
        multiqc -n QC_report.html .
        """
}
