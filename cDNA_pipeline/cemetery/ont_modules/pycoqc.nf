process PYCOQC {

    publishDir 'results/pycoQC/', mode: 'copy', overwrite: false

    label 'large'

    input:
        val(id)
        path(seq_summary)
        path(total_bam)
        path(total_bai)

    output:
        path "${id}_pycoqc.json", emit: multiqc
        path "*", emit: all_output
    
    script:
    """
    fix_pyco_summary.py $seq_summary "final_pyco.txt"

    pycoQC -f "final_pyco.txt" \
        -a $total_bam \
        -o "./${id}_pycoqc.html" \
        -j "./${id}_pycoqc.json" \
    """ 
}

