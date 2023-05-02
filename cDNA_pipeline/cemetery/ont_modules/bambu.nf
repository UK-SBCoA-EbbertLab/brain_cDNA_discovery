process BAMBU {

    publishDir 'results/bambu/', mode: 'copy', overwrite: false

    label 'large'

    input:
        val(id)
        tuple val(d), path(bam)
        tuple val(d2), path(bai)
        path(ref)
        path(gtf)
        path(fai)
        

    output:
        path("*")

    shell:
    '''
    mkdir merged_final/
    mkdir merged_rds/

    bam_dummy="!{bam}"

    bam2="$(tr ' ' ',' <<<$bam_dummy)"
    
    echo $bam2 > test_file.txt

    bambu_script.R $bam2 "!{ref}" "!{gtf}" "!{id}"
    '''
}

// POSSIBLY TAKE ID VAR WAY FROM BAMBU SCRIPT
