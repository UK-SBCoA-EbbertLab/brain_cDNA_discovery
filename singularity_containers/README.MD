# Singularity Containers




bambu.def - definition file for singularity container used to run the Bambu R package on the Nextflow Pipeline.

`pull command: singularity pull --arch amd64 library://bernardo-a-heberle/cdna_analysis/bambu:sha256.0e14bc0ebf3ca9373d0258276776710dffc156cbb7c62b0eb891157c7d90254b`



bernardo_article_analysis.def - definition file for singularity container used to perform the analysis downstream from the NextFlow Pipeline. 

`pull command: library://bernardo-a-heberle/cdna_analysis/bernardo_article_analysis:sha256.7b9f63b038fbf11c4b5aa68eacc96a5058a2a20d74f2a239f0254b335bf4ffda`



guppy.def - definition file for singularity container used to run guppy basecaller... Was not used due to data being basecalled on the PromethION.

 `pull command: library://bernardo-a-heberle/cdna_analysis/guppy:sha256.9a9819832e237723a91274e81a8c803bde20b46e8393c610ad71555cf148065f`
 
 
 
meme_suite.def - definition file for singularity container used to run the meme-suite package.

 `pull commmand: library://bernardo-a-heberle/cdna_analysis/meme_suite:sha256.cd6370655c228734af335070a6c5195627fea6c6314e1ca066a66f2584b992da`
 



nanopore.def - definition file for singularity container used to run the software nanopore data analysis.

`pull command: library://bernardo-a-heberle/cdna_analysis/nanopore:sha256.c37569b9b0c36cb71d571a7c5ac8580a8aee439e0cea4a7a54e754dad0191225`



quality_control.def - definition file for singularity container used to run the quality control software.

`pull command: singularity pull --arch amd64 library://bernardo-a-heberle/cdna_analysis/quality_control:sha256.527a5b5f4dab01c393b08456a0f0188378d9603fd253f587d63a6171a3bc421d`




### For more information about software versions see the %help section of definition (.def) files for the singularity images.
