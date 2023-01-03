#!/bin/bash
#SBATCH -p defq
#SBATCH -n 3
#SBATCH --time=25:00:00
#SBATCH -o 9.caseout.txt
#SBATCH -e 9.caseerr.txt
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

matlab -r build_mutation_table_master_fast
echo Done!!!