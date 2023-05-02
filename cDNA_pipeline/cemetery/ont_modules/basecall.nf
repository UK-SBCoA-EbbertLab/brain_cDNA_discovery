process BASECALL_RNA {

    label 'gpu'

    input:
    tuple val(read_num), path(fast5_dir)

    output:
    path('pass/*.fastq'), emit: fastq
    path('*.txt'), emit: txt

   script:
    """
    guppy_basecaller --device cuda:all:100% --input_path '${fast5_dir}' -q 100000000000000 --save_path './' -c rna_r9.4.1_70bps_hac_prom.cfg
    
    mv pass/*.fastq 'pass/${read_num}.fastq'
    mv *.txt '${read_num}.txt'
    """
}

process GATHER_BASECALL {

    publishDir 'results/guppy/', mode: "copy", overwrite: true

    label 'tiny'
 
    input:
    val(id)
    path('*.fastq')
    path('*.txt')

    output:
    path('*.fastq'), emit: fastq
    path('*.txt'), emit: txt

    script:
    """
    cat *.fastq > "${id}.fastq"
    cat *.txt > "${id}.txt"
    """
}
