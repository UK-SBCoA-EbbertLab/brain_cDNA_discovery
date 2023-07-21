#!/bin/bash


## Start jupyter notebooks through singularity container and tunnel it through port 1234
## Need to have used ssh <user@serve> -L 1234:localhost:1234 for the tunneling to work.
singularity exec ../../singularity_containers/bernardo_article_analysis.sif jupyter notebook --no-browser --port 1234
