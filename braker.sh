#!/bin/bash
#SBATCH -A e30514
#SBATCH -p normal
#SBATCH --job-name="snowfly braker"
#SBATCH -N 1
#SBATCH --ntasks-per-node=52
#SBATCH -t 48:00:00
#SBATCH --mem=180G

cd $SLURM_SUBMIT_DIR

module purge all
module load singularity
export BRAKER_SIF=/home/rsd3363/software/braker3.sif

# Here, we run the BRAKER3 docker image through Singularity for ease of use 
# with a cluster. Inputs are a protein library under "--prot_seq", read files 
# bam format under "--bam", and a config path for augustus under "--AUGUSTUS_CONFIG_PATH"


# OrthoDB clades (used for protein library) were partitioned by using this code
# https://github.com/tomasbruna/orthodb-clades


# In case you rerun the script since if this file exists, BRAKER spits out error
rm -r /home/rsd3363/.augustus/species/snowfly 


# Run BRAKER3 with singularity. Use 4 less threads in the command than you allocate. 
singularity exec -B /projects/b1020/Richard:/Richard ${BRAKER_SIF} braker.pl --species=snowfly --genome=/Richard/snowfly/Chionea_cleaned_v2.softmasked.fasta --prot_seq=/Richard/other_data/proteins/orthodb11/orthodb-clades/clades/Diptera.fa \
--bam=/Richard/snowfly/align/snowfly_ab_f_Aligned.sortedByCoord.out.bam,/Richard/snowfly/align/snowfly_ab_m_Aligned.sortedByCoord.out.bam,/Richard/snowfly/align/snowfly_ant_Aligned.sortedByCoord.out.bam,/Richard/snowfly/align/snowfly_head_Aligned.sortedByCoord.out.bam,/Richard/snowfly/align/snowfly_thorax_Aligned.sortedByCoord.out.bam \
--threads 48 --AUGUSTUS_CONFIG_PATH=/home/rsd3363/.augustus --workingdir=/Richard/snowfly/braker


# Convert gtf to gff for Mikado using AGAT suite
agat_convert_sp_gxf2gxf.pl -g braker.gtf -o braker.gff



