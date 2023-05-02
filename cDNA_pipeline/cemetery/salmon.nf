process SALMON_ALIGNMENT_MODE {

    publishDir "results/${params.out_dir}/salmon_alignment_mode/"

    label "medium"

    input:
        val(id)
        val(d)
        path(bam)
        path(transcriptome)

    output:
        path "*" 


    script:
        """ 
        salmon quant -t $transcriptome -l A -a $bam -o '${id}_${d}_salmon_alignment_based_quant'
        """ 

}
