#!/bin/bash
#SBATCH -p defq,sched_mem1TB
#SBATCH -n 16
#SBATCH --time=12:00:00
#SBATCH -o 7.masterout.txt
#SBATCH -e 7.mastererr.txt
#SBATCH --mem=80000

module add c3ddb/bamtools/2.4.0
module add c3ddb/bcftools/1.2
module add c3ddb/bowtie2/2.2.6
module add c3ddb/htslib/1.2.1
module add c3ddb/samtools/1.2
module add c3ddb/sickle/1.33
module add c3ddb/samtools/1.2
module add mit/matlab/2015b
module add c3ddb/dindel-tgi

matlab -r build_experiment_directories_unpaired
echo Done!!!
