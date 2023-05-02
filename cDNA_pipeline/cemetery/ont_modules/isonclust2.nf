process ISONCLUST2 {
    
    conda 'isONclust2/env.yml'
    publishDir 'results/isONclust2/', mode: 'copy', overwrite: false
    
    label "big_mem"

    input:
        val(id)
        path fastq

    output:
        path "output_directory"   

    script:
    """ 
    . /conda/etc/profile.d/conda.sh
    conda activate denovo-isoforms
    git clone --recursive https://github.com/bernardo-heberle/pipeline-nanopore-denovo-isoforms
    cd pipeline-nanopore-denovo-isoforms
    snakemake -j 4 all
    """ 
}

