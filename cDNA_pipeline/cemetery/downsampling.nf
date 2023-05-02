process DOWNSAMPLING_ILLUMINA {

    publishDir "results/${params.out_dir}/illumina_downsampling/", mode: "copy", overwrite: true

    label 'medium_small'

    input:
        val(id)
        path(bam_trans)
        path(bam_gen)
        each d

    output:
        val "$id", emit: id
        path "${id}_trans_downsampled_${d}.bam", emit: bam_salmon
        path "${id}_gen_downsampled_${d}_sorted.bam", emit: bam_stringtie
        path "${id}_gen_downsampled_${d}_sorted.bam.bai", emit: bai_stringtie
        val "$d", emit: downsample

    script:
        """
        samtools view -@ 12 -s $d -b $bam_trans > "${id}_trans_downsampled_${d}.bam"
        
        samtools view -@ 12 -s $d -b $bam_gen > "${id}_gen_downsampled_${d}.bam"
        samtools sort -@ 12 "${id}_gen_downsampled_${d}.bam" -o "${id}_gen_downsampled_${d}_sorted.bam"
        samtools index -@ 12 "${id}_gen_downsampled_${d}_sorted.bam"
        """
}



process DOWNSAMPLING_NANOPORE {

    publishDir "results/${params.out_dir}/nanopore_downsampling/", mode: "copy", overwrite: true

    label 'medium_small'

    input:
        val(id)
        path(bam)
        each d

    output:
        val "$id", emit: id
        path "${id}_downsampled_${d}.bam", emit: bam
        path "${id}_downsampled_${d}.bam.bai", emit: bai
        val "$d", emit: downsample

    script:
        """
        samtools view -@ 12 -s $d -b $bam > "${id}_downsampled_${d}.bam"
        samtools index -@ 12 "${id}_downsampled_${d}.bam"
        """
}
