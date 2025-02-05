#!/bin/bash
#SBATCH -A b1020
#SBATCH -p b1020
#SBATCH --job-name="snowfly align"
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH -t 24:00:00
#SBATCH --mem=80G

module load STAR/2.7.9a
module load samtools
cd $SLURM_SUBMIT_DIR

# Here we align the data to STAR. In theory, this might need to be done twice. 
# The first time will be for Blobtools to look at coverage. If the assembly is
# sufficiently clean, then you can move on. Otherwise, remove "bad" scaffolds
# and rerun STAR. Here, we run STAR with the options 
# "--outSAMstrandField intronMotif" and "--outFilterIntronMotifs RemoveNoncanonical" 
# as they are reccomended for running Stringtie. The other options are specific 
# to species/read file type, etc. 


MAXINTRON=50000

DNA_LOC=/projects/b1020/Richard/snowfly/data/Chionea_cleaned_v2.softmasked.fasta
GENOME_DIR=/projects/b1020/Richard/snowfly/indices/index
RF_PREFIX=/projects/b1020/Richard/snowfly/fastq/

# Build STAR database, this only needs to be done once so if the following 
# STAR commands fail, comment this one out
STAR --runMode genomeGenerate --runThreadN 20 --genomeDir $GENOME_DIR --genomeFastaFiles $DNA_LOC --genomeSAindexNbases 13

# Run STAR for each fastq file (or each pair in some cases)

# ab_f
STAR --genomeDir $GENOME_DIR --readFilesPrefix $RF_PREFIX --readFilesIn ABDOMEN_F_S5_R1_001.fastq.gz --readFilesCommand "gzip -cd" --outFileNamePrefix snowfly_ab_f_ --runThreadN $SLURM_NTASKS_PER_NODE --outFilterType BySJout --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignIntronMax $MAXINTRON --alignMatesGapMax $MAXINTRON --outFilterIntronMotifs RemoveNoncanonical 

samtools index snowfly_ab_f_Aligned.sortedByCoord.out.bam

# ab_m
STAR --genomeDir $GENOME_DIR --readFilesPrefix $RF_PREFIX --readFilesIn ABDOMEN_M_S4_R1_001.fastq.gz --readFilesCommand "gzip -cd" --outFileNamePrefix snowfly_ab_m_ --runThreadN $SLURM_NTASKS_PER_NODE --outFilterType BySJout --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignIntronMax $MAXINTRON --alignMatesGapMax $MAXINTRON --outFilterIntronMotifs RemoveNoncanonical 

samtools index snowfly_ab_m_Aligned.sortedByCoord.out.bam

# ant
STAR --genomeDir $GENOME_DIR --readFilesPrefix $RF_PREFIX --readFilesIn ANTENNAE_S2_R1_001.fastq.gz --readFilesCommand "gzip -cd" --outFileNamePrefix snowfly_ant_ --runThreadN $SLURM_NTASKS_PER_NODE --outFilterType BySJout --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignIntronMax $MAXINTRON --alignMatesGapMax $MAXINTRON --outFilterIntronMotifs RemoveNoncanonical 

samtools index snowfly_ant_Aligned.sortedByCoord.out.bam

# head
STAR --genomeDir $GENOME_DIR --readFilesPrefix $RF_PREFIX --readFilesIn HEADS_S1_R1_001.fastq.gz --readFilesCommand "gzip -cd" --outFileNamePrefix snowfly_head_ --runThreadN $SLURM_NTASKS_PER_NODE --outFilterType BySJout --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignIntronMax $MAXINTRON --alignMatesGapMax $MAXINTRON --outFilterIntronMotifs RemoveNoncanonical 

samtools index snowfly_head_Aligned.sortedByCoord.out.bam

# thorax
STAR --genomeDir $GENOME_DIR --readFilesPrefix $RF_PREFIX --readFilesIn THORAX_S3_R1_001.fastq.gz --readFilesCommand "gzip -cd" --outFileNamePrefix snowfly_thorax_ --runThreadN $SLURM_NTASKS_PER_NODE --outFilterType BySJout --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignIntronMax $MAXINTRON --alignMatesGapMax $MAXINTRON --outFilterIntronMotifs RemoveNoncanonical 

samtools index snowfly_thorax_Aligned.sortedByCoord.out.bam
