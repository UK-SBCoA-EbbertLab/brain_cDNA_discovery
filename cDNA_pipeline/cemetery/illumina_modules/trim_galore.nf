process DECOMPRESS {

    label "tiny"

    input:
    tuple val(id), path(file_R1), path(file_R2)

    output:
    val "$id", emit: id
    path "*_R1*.fastq", emit: file_R1
    path "*_R2*.fastq", emit: file_R2

    script:
    """
    gzip -f -d $file_R1

    gzip -f -d $file_R2
    """

}

process TRIM_GALORE {

    publishDir 'results/trim_galore/', mode: 'copy', overwrite: false

    label "small"

    input:
    val(id)
    path(file_R1)
    path(file_R2)

    output:
    val "$id", emit: id
    path "*val_1*.fq", emit: trim_1
    path "*val_2*.fq", emit: trim_2
    path "*report.txt", emit: QC_1
    path "*fastqc.zip", emit: QC_2

    script:
    """
    trim_galore --paired --illumina -j 8 --fastqc "${file_R1}" "${file_R2}" -o ./ --basename "${id}"
    """
}
