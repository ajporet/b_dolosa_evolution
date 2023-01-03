#!/bin/bash
#SBATCH -p defq,sched_mem1TB
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=72:00:00
#SBATCH -o masteroutA.txt
#SBATCH -e mastererrA.txt
#SBATCH --job-name="1_get_SRA_files"
#SBATCH --mem=8000
#SBATCH --array=1-58



module add c3ddb/sratoolkit/2.8.0

sraID=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /scratch/mit_lieberman/projects/gga_bdelosa/09_22_Past_2011_Data/sra_accessions_partA.txt |cut -f1)


fastq-dump -O /scratch/mit_lieberman/projects/gga_bdelosa/09_22_Past_2011_Data/fastq_files/ -I --gzip ${sraID}