#!/bin/bash
#SBATCH -p defq,sched_mem1TB
#SBATCH -n 1
#SBATCH --time=30:00:00
#SBATCH -o masterout.txt
#SBATCH -e mastererr.txt
#SBATCH --mem=8000
#SBATCH --job-name="breseq_snakemake"
bash snakemakeslurm.sh
echo Done!!!
