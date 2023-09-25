process MULTIQC_GRCh38 {

    publishDir "results/${params.out_dir}/QC/multiqc", mode: "copy", overwrite: true

    label 'tiny'

    input:
        path(QC_1)
        path(QC_2)
        path(QC_3)
        path(QC_4)
        path(QC_5)
        path(QC_6)
        path(QC_7)
        path(multiqc_config)
    
    output: 
        path "*"

    script:
        """    
        multiqc -c $multiqc_config -n multiQC_report.html .
        """
}

process MULTIQC_CHM13 {

    publishDir "results/${params.out_dir}/QC/multiqc/", mode: "copy", overwrite: true

    label 'tiny'

    input:
        path(QC_1)
        path(QC_2)
        path(QC_3)
        path(QC_4)
        path(QC_5)
        path(QC_6)
        path(multiqc_config)
    
    output: 
        path "*"

    script:
        """    
        multiqc -c $multiqc_config -n multiQC_report.html .
        """
}

process MULTIQC_GRCh38_ONT_ONLY {

    publishDir "results/${params.out_dir}/QC/multiqc", mode: "copy", overwrite: true

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
        multiqc -c $multiqc_config -n multiQC_report.html .
        echo "hi"
        """
}

process MULTIQC_CHM13_ONT_ONLY {

    publishDir "results/${params.out_dir}/QC/multiqc", mode: "copy", overwrite: true

    label 'tiny'

    input:
        path(QC_1)
        path(QC_2)
        path(QC_3)
        path(multiqc_config)
    
    output: 
        path "*"

    script:
        """    
        multiqc -c $multiqc_config -n multiQC_report.html .
        """
}
