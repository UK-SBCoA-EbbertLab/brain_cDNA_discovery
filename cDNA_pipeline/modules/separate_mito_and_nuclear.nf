process SEPARATE_MITO_ANNOTATION {

    publishDir "results/${params.out_dir}/separated_annotations/", mode: "copy", overwrite: true

    label 'medium'

    input:
        path fai
        path fasta
        path gtf

    output:
        path('nuclear_ref.fa'), emit: nuclear_fa
        path('mito_ref.fa'), emit: mito_fa
        path('nuclear_ref.fa.fai'), emit: nuclear_fai
        path('mito_ref.fa.fai'), emit: mito_fai
        path('nuclear_annotation.gtf'), emit: mito_gtf
        path('mito_annotation.gtf'), emit: nuclear_gtf

    shell:
        '''
        for chr in $(cut -f1 "!{fai}")
        do
            if [ "$chr" = "chrM" ] || [ "$chr" = "MT" ]; then
                samtools faidx "!{fasta}" $chr > "mito_ref.fa"
                mito_chr_name=$chr
            else
                samtools faidx "!{fasta}" $chr >> "nuclear_ref.fa"
            fi
        done

        samtools faidx nuclear_ref.fa
        samtools faidx mito_ref.fa
        
        awk -F '\t' '$1=="${mito_chr_name}"' "!{gtf}" >> mito_annotation.gtf
        awk -F '\t' '$1!="${mito_chr_name}"' "!{gtf}" >> nuclear_annotation.gtf

        '''

}
