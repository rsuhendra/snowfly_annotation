#!/bin/bash
#SBATCH -A p31669
#SBATCH -p normal
#SBATCH --job-name="snowfly repeatmodel"
#SBATCH -N 1
#SBATCH --ntasks-per-node=52
#SBATCH -t 48:00:00
#SBATCH --mem=80G
#SBATCH --output=repeatmodeler.log
#SBATCH --mail-user=rsuhendra@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR

module purge all
source activate augusco

# Here we run RepeatModeler to build up a de-novo library of repeats to 
# be used for RepeatMasker

# Build RepeatModeler Database
BuildDatabase -name snowfly -engine ncbi /projects/b1020/Richard/snowfly/data/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.fasta

# Run RepeatModeler
RepeatModeler -database snowfly -threads 52 -LTRStruct

# Adds the prefix "chiAle1" to each fasta header
cat snowfly-families.fa | seqkit fx2tab | awk '{ print "chiAle1_"$0 }' | seqkit tab2fx > snowfly-families.prefix.fa

# seperate known and unknown elements
cat snowfly-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > snowfly-families.prefix.fa.known

cat snowfly-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > snowfly-families.prefix.fa.unknown