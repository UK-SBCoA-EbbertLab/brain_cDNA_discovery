process DOWNSAMPLING {

    publishDir 'results/downsampling/', mode: 'copy', overwrite: false

    label 'medium_small'

    input:
        val(id)
        path(bam)
        each d

    output:
        val "$id", emit: id
        val "${d}", emit: downsample
        path "${id}_downsampled_${d}.bam", emit: bam 
    
    script:
    """
    samtools view -@ 12 -s $d -b $bam > "${id}_downsampled_${d}.bam"
    """
}
