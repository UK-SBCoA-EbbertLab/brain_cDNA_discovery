process TRIM_GALORE {

    publishDir "results/${params.out_dir}/trim_galore/"
    publishDir "results/${params.out_dir}/QC/fastqc/", mode: 'copy', overwrite: true, pattern: "*fastqc.zip*"

    label "small"

    input:
        tuple val(id), path(file_R1), path(file_R2)

    output:
        val "$id", emit: id
        path "*val_1*.fq", emit: trim_1
        path "*val_2*.fq", emit: trim_2
        path "*report.txt", emit: QC
        path "*fastqc.zip", emit: multiQC

    script:
        """
        trim_galore --paired --dont_gzip --illumina -j 8 --fastqc "${file_R1}" "${file_R2}" -o ./ --basename "${id}"
        """
}
