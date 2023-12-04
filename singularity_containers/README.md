# Singularity Containers




bambu.def - definition file for singularity container used to run the Bambu R package on the Nextflow Pipeline.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/bambu:sha256.c766826dd183d2c09be2ae4b64524954cecc66ea31557483fda249dd34c21c1d`



bernardo_article_analysis.def - definition file for singularity container used to perform the analysis downstream from the NextFlow Pipeline. 

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/bernardo_article_analysis:sha256.c44e2b51cd227900f9e3df917dcef07a70daf0339baf58d26dce3b638e54b191`



guppy.def - definition file for singularity container used to run guppy basecaller... Was not used due to data being basecalled on the PromethION.

 `pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/guppy:sha256.80d73ed421e17199390658662453c2c37c8e64435d8b402822a107247882915f`


 illumina.def - definition file for singularity container used to run the `illumina_pipeline` NextFlow pipeline.
 `pull command: singularity pull --arch amd64 library://bernardo-a-heberle/cdna_analysis/illumina:sha256.271989b636e157fe3032a1339fb3d632b6efe88e8726e79e3d76cdd038035f77 `
 
 
 
meme_suite.def - definition file for singularity container used to run the meme-suite package.

 `pull commmand: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/meme_suite:sha256.8696beb153ef66f4966aa9ee5528cf7b0b4da1f10f2915a67be44cbc225cc781`
 



nanopore.def - definition file for singularity container used to run the software nanopore data analysis.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/nanopore:sha256.e8cb0f02b0e8b3f89937e88333f3c410eff3eacc53d861f5c04a1768790c6a21`



quality_control.def - definition file for singularity container used to run the quality control software.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/quality_control:sha256.2db86def19a30f36a55ffbf8f999803cace938388f13a0c861c3f75efe49f1e7`




### For more information about software versions see the %help section of definition (.def) files in this directory.
