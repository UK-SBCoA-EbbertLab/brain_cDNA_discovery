process BASECALL {

    label 'gpu'

    input:
        path fast5
        val config

    output:
        path 'pass/*.fastq', emit: fastq
        path '*.txt', emit: txt

   script:
        """
        guppy_basecaller --input_path '!{fast5}' --save_path './' -c $config
        """
}

process GATHER_BASECALL {

    publishDir "results/${params.out_dir}/basecall_output/", mode: "copy", overwrite: true

    label 'local'

    input:
        val id
        path fastq
        path txt

    output:
        path 'total_basecall.fastq', emit: fastq
        path 'total_basecall.txt', emit: txt
        path '*', emit: output

    shell:
        '''
        cat "!{fastq}" >> "!{id}_total_basecall.fastq"
        cat "!{txt}" >> "!{id}_total_basecall.txt"
        '''
}
