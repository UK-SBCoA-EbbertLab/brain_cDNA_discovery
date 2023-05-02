process DOWNSAMPLING {

    publishDir 'results/downsampling/', mode: 'copy', overwrite: false

    label 'medium_small'

    input:
        val(id)
        path(bam)
        each d

    output:
        val "$id", emit: id
        tuple val("$d"), path("${id}_downsampled_${d}.bam"), emit: bam
        tuple val ("$d"), path("${id}_downsampled_${d}.bam.bai"), emit: bai

    script:
    """
    samtools view -@ 12 -s $d -b $bam > "${id}_downsampled_${d}.bam"
    samtools index -@ 12 "${id}_downsampled_${d}.bam"
    """
}
