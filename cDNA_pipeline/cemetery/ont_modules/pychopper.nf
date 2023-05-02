process PYCHOPPER {

    publishDir 'results/pychopper/', mode: 'copy', overwrite: false

    label "large"

    input:
        tuple val(id), path(fastq)

    output:
        val "$id", emit: id
        path "${id}_pychop.fq", emit: fastq
        path "*pychopper*", emit: multiqc

    script:
    """
    cdna_classifier.py -t 50 \
        -r "${id}_pychopper_report.pdf" \
        -u "${id}_pychopper.unclassified.fq" \
        -w "${id}_pychopper.rescued.fq" \
        -S "${id}_pychopper.stats" \
        -A "${id}_pychopper.scores" \
        "${fastq}" "${id}_pychop.fq"
    """
}
