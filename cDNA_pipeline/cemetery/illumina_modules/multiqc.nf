process MULTIQC {

    publishDir 'results/multiqc/', mode: 'copy', overwrite: false

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
    
    multiqc -n dcnl_illumina_QC_report.html .

    """



}
