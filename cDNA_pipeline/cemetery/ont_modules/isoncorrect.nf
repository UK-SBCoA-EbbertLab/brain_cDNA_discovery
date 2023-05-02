process ISONCORRECT {
    
    publishDir 'results/isONcorrect/', mode: 'copy', overwrite: false
    
    label "pychopper"

    input:
        val(id)
        path fastq

    output:
        val "$id", emit: id
        path "*all_corrected_reads.fq", emit: fastq   

    script:
    """
    bash isoncorrect.sh $fastq $id 62
    """
}


process ISONCORRECT_PARALELIZED_RNA {

    publishDir 'results/isONcorrect/', mode: 'copy', overwrite: false
    
    label "pychopper"

    input:
        val(id)
        path(fastq)
        each i

    output:
       tuple val("$id"), path("*corrected_reads.fq"), emit: fastq
       path "*", emit: output

    script:
    """
    isoncorrect_paralelized.sh $fastq $id 62 $i
    """
}


process ISONCORRECT_PARALELIZED_CDNA {

    publishDir 'results/isONcorrect/', mode: 'copy', overwrite: false

    label "pychopper"

    input:
        val(id)
        path(fastq)
        each i

    output:
       tuple val("$id"), path("*corrected_reads.fq"), emit: fastq
       path "*", emit: output

    script:
    """
    isoncorrect_paralelized.sh $fastq $id 62 $i
    """
}


process ISONCLUST {

    publishDir 'results/isonclust/', mode: 'copy', overwrite: false

    label "pychopper"

    input:
        val(id)
        path(fastq)

    output:
        val "$id", emit: id
        path "clustering/fastq_files/", emit: fastq

    script:
    """
    bash isonclust.sh $fastq $id 50
    """
}

process ISONGATHER{
    
    publishDir 'results/corrected_reads/', mode: 'copy', overwrite: false
    
    label "tiny"
    
    input:
        tuple val(id), path(fastq)
       
    output:
        val "$id", emit: id
        path "${id}_corrected.fq", emit: fastq

    script:
    """
    cat $fastq >> "${id}_corrected.fq"
    """
}


process ISONCLUST_RNA {


    publishDir 'results/isonclust/', mode: 'copy', overwrite: false

    label "pychopper"

    input:
        tuple val(id), path(fastq)

    output:
        val "$id", emit: id
        path "clustering/fastq_files/", emit: fastq

    script:
    """
    bash isonclust.sh $fastq $id 50
    """
}
