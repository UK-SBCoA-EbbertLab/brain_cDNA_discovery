#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --job-name=MSV000086944_fragpipe
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu 25G
#SBATCH --partition=normal
#SBATCH -e slurm-%j.err           
#SBATCH -o slurm-%j.out            
#SBATCH -A coa_mteb223_uksr

module load ccs/java/jdk-19.0.2
cd /scratch/bjwh228/fragpipe/bin
./fragpipe --headless --workflow /scratch/bjwh228/fragpipe/workflows/fragpipe.workflow --manifest /scratch/bjwh228/fragpipe/MSV000086944_noFail_manifest.fp-manifest --workdir /scratch/bjwh228/results_FDR_10_perc_ALL_for_cell_lines/ --config-msfragger /scratch/bjwh228/fragpipe/tools/MSFragger-3.8/MSFragger-3.8.jar --config-ionquant /scratch/bjwh228/fragpipe/tools/IonQuant-1.9.8.jar --config-philosopher /scratch/bjwh228/fragpipe/tools/philosopher_v5.0.0_linux_amd64/philosopher --config-python /share/apps/amd/conda/envs/python-3.9.6/bin/python3.9
