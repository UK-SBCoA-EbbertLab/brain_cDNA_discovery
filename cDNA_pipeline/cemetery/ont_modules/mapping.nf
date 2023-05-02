process MAPPING {

    publishDir 'results/mapping/', mode: 'copy', overwrite: false

    label 'large'

    input:
        val(id)
        path(fastq)
        path(index)

    output:
        val "$id", emit: id
        path "${id}_mapped_sorted.bam", emit: bam
        path "${id}_mapped_sorted.bam.bai", emit: bai
        path "*stat", emit: multiqc
        path "*", emit: output

    script:
    """
    echo id
    minimap2 -t 16 -ax splice \
        -uf \
        $index \
        $fastq > "${id}_all.bam" \
     
    samtools view -b -F 4 -@ 12 '${id}_all.bam' > '${id}_mapped.bam'
    samtools sort -@ 16 "${id}_mapped.bam" -o '${id}_mapped_sorted.bam'
    samtools index '${id}_mapped_sorted.bam'
    samtools flagstat "${id}_mapped_sorted.bam" > "${id}_mapped_sorted.flagstat"
    samtools idxstats "${id}_mapped_sorted.bam" > "${id}_mapped_sorted.idxstat"
     """
}

process MAPPING_QC {

    publishDir 'results/mapping_no_correction/', mode: 'copy', overwrite: false

    label 'large'

    input:
        tuple val(id), path(fastq)
        path(index)

    output:
        val "$id", emit: id
        path "${id}_mapped_sorted.bam", emit: bam
        path "${id}_mapped_sorted.bam.bai", emit: bai
        path "*", emit: output

    script:
    """ 
    minimap2 -t 16 -ax splice \
        -uf \
        $index \
        $fastq > '${id}_all.bam' \
     
    samtools view -b -F 4 -@ 12 '${id}_all.bam' > '${id}_mapped.bam'
    samtools sort -@ 16 "${id}_mapped.bam" -o '${id}_mapped_sorted.bam'
    samtools index '${id}_mapped_sorted.bam'
    samtools flagstat "${id}_mapped_sorted.bam" > "${id}_mapped_sorted.flagstat"
    samtools idxstats "${id}_mapped_sorted.bam" > "${id}_mapped_sorted.idxstats"
     """ 
}

process MAPPING_QC_DIRECT_RNA {

    publishDir 'results/mapping_no_correction_direct/', mode: 'copy', overwrite: false

    label 'large'

    input:
        tuple val(id), path(fastq)
        path(index)

    output:
        val "$id", emit: id
        path "${id}_all_sorted.bam", emit: bam_all
        path "${id}_all_sorted.bam.bai", emit: bai_all
        path "${id}_mapped_sorted.bam", emit: bam_mapped
        path "${id}_mapped_sorted.bam.bai", emit: bai_mapped
        path "*stat", emit: multiqc
        path "*", emit: output


    script:
    """
    minimap2 -t 16 -ax splice \
        -k14 \
        $index \
        $fastq > '${id}_all.bam' \
    
     
    samtools sort -@ -16 "${id}_all.bam" -o "${id}_all_sorted.bam"
    samtools index "${id}_all_sorted.bam"
    samtools flagstat "${id}_all_sorted.bam" > "${id}_all_sorted.flagstat"
    samtools idxstats "${id}_all_sorted.bam" > "${id}_all_sorted.idxstat"

    samtools view -b -F 4 -@ 12 '${id}_all.bam' > '${id}_mapped.bam'
    samtools sort -@ 16 "${id}_mapped.bam" -o '${id}_mapped_sorted.bam'
    samtools index '${id}_mapped_sorted.bam'
     """
}

process MAPPING_DIRECT_RNA {

    publishDir 'results/mapping_correction_direct/', mode: 'copy', overwrite: false

    label 'large'

    input:
        val(id)
        path(fastq)
        path(index)

    output:
        val "$id", emit: id
        path "${id}_all_sorted.bam", emit: bam_all
        path "${id}_all_sorted.bam.bai", emit: bai_all
        path "${id}_mapped_sorted.bam", emit: bam_mapped
        path "${id}_mapped_sorted.bam.bai", emit: bai_mapped
        path "*stat", emit: multiqc
        path "*", emit: output

    script:
    """
    minimap2 -t 16 -ax splice \
        -k14 \
        $index \
        $fastq > '${id}_all.bam' \

    
    samtools sort -@ -16 "${id}_all.bam" -o "${id}_all_sorted.bam"
    samtools index "${id}_all_sorted.bam"
    samtools flagstat "${id}_all_sorted.bam" > "${id}_all_sorted.flagstat"
    samtools idxstats "${id}_all_sorted.bam" > "${id}_all_sorted.idxstat"

    samtools view -b -F 4 -@ 12 '${id}_all.bam' > '${id}_mapped.bam'
    samtools sort -@ 16 "${id}_mapped.bam" -o '${id}_mapped_sorted.bam'
    samtools index '${id}_mapped_sorted.bam'
    """

}
