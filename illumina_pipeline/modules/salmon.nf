process SALMON_ALIGNMENT_MODE {

    publishDir "results/${params.out_dir}/salmon_alignment_mode/", mode: "copy", overwrite: true

    label "medium"

    input:
        val(id)
        path(bam)
        path(transcriptome)

    output:
        path "*" 


    script:
        """ 
        salmon quant -t $transcriptome -l A -a $bam -o '${id}_salmon_alignment_based_quant'
        """ 

}
