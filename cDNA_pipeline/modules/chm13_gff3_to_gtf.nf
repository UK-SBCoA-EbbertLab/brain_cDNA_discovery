process CHM13_GTF_ERCC {

    publishDir "results/${params.out_dir}/CHM13_gtf/", mode: "copy", overwrite: true

    label 'medium'

    input:
        path gff
        path ercc

    output:
        path('CHM13_v2.0_ERCC.gtf')

    script:
        """
        gff_to_gtf.py $gff CHM13_v2.0.gtf

        cat CHM13_v2.0.gtf $ercc > CHM13_v2.0_ERCC.gtf
        """
}

process CHM13_GTF {

    publishDir "results/${params.out_dir}/CHM13_gtf/", mode: "copy", overwrite: true

    label 'medium'

    input:
        path gff

    output:
        path('CHM13_v2.0.gtf')

    script:
        """
        gff_to_gtf.py $gff CHM13_v2.0.gtf
        """
}

