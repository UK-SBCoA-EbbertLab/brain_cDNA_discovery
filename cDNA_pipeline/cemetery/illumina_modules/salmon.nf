process SALMON_MAPPING_MODE {

    publishDir 'results/salmon_mapping_mode/', mode: 'copy', overwrite: false

    label "medium"

    input:
    val(id)
    path(file_R1)
    path(file_R2)
    path(gentrome)
    path(decoys)

    output:
    path "*"


    script:
    """
    salmon index -t $gentrome -i transcripts_index -d $decoys -k 31 -p 12 
    
    salmon quant -i transcripts_index -l A -1 $file_R1 -2 $file_R2 --validateMappings -o ${id}_salmon_mapping_based_quant
    """

}

process SALMON_ALIGNMENT_MODE {


    publishDir 'results/salmon_alignment_mode/', mode: 'copy', overwrite: false

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
