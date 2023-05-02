process FIX_SEQ_SUMMARY {

    label "medium"

    input:
        val(id)
        path(seq_summary)
        path(fastq_pychop)

    output:
        val "$id", emit: id
        path "${id}_sequencing_summary_fixed.txt", emit: txt
        path "$fastq_pychop", emit: fastq


    script:
        """
        fix_sequencing_summary_pychopper.py $fastq_pychop $seq_summary "${id}_sequencing_summary_fixed.txt"
        """
}
