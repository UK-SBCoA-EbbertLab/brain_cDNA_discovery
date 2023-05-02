process PYCHOPPER {

    publishDir "results/${params.out_dir}/pychopper/"
    publishDir "results/${params.out_dir}/QC/pychopper/", mode: 'copy', overwrite: true, pattern: "*pychopper.stats"

    label "large"

    input:
        tuple val(id), path(fastq)
        path(txt)
        val(cdna_kit)

    output:
        val "$id", emit: id
        path "${id}_pychop.fq", emit: fastq
        path "$txt", emit: txt
        path "$fastq", emit: original_fastq
        path "*pychopper.stats", emit: multiQC

    script:
        """
        pychopper -t 50 \
            -Q 9 \
            -k $cdna_kit \
            -r "${id}_pychopper_report.pdf" \
            -u "${id}_pychopper.unclassified.fq" \
            -w "${id}_pychopper.rescued.fq" \
            -S "${id}_pychopper.stats" \
            -A "${id}_pychopper.scores" \
            "${fastq}" "${id}_pychop.fq"
        """
}
