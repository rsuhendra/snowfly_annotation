#!/bin/bash
#SBATCH -A b1020
#SBATCH -p b1020
#SBATCH --job-name="stringtie snowfly"
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -t 4:00:00
#SBATCH --mem=40G

module load stringtie

cd $SLURM_SUBMIT_DIR

# Here, we use run Stringtie, Strawberry, and Scallop, which are transcriptome generating
# tools using only aligned RNA seq reads. We run each program once for each bam file. 
# We use mostly default settings except in the case of Stringtie, where we elected to 
# use the "--conservative" option. 

# Stringtie commands, multi-thread
stringtie /projects/b1020/Richard/snowfly/alignments/align/snowfly_ab_f_Aligned.sortedByCoord.out.bam --conservative -o snowfly_stringtie_ab_f.gtf -p 12 -A snowfly_stringtie_ab_f_abundance.tab

stringtie /projects/b1020/Richard/snowfly/alignments/align/snowfly_ab_m_Aligned.sortedByCoord.out.bam --conservative -o snowfly_stringtie_ab_m.gtf -p 12 -A snowfly_stringtie_ab_m_abundance.tab

stringtie /projects/b1020/Richard/snowfly/alignments/align/snowfly_head_Aligned.sortedByCoord.out.bam --conservative -o snowfly_stringtie_head.gtf -p 12 -A snowfly_stringtie_head_abundance.tab

stringtie /projects/b1020/Richard/snowfly/alignments/align/snowfly_ant_Aligned.sortedByCoord.out.bam --conservative -o snowfly_stringtie_ant.gtf -p 12 -A snowfly_stringtie_ant_abundance.tab

stringtie /projects/b1020/Richard/snowfly/alignments/align/snowfly_thorax_Aligned.sortedByCoord.out.bam --conservative -o snowfly_stringtie_thorax.gtf -p 12 -A snowfly_stringtie_thorax_abundance.tab

# Strawberry commands, multi-thread
strawberry ../alignments/align/snowfly_ab_f_Aligned.sortedByCoord.out.bam -o snowfly_ab_f_strawberry.gtf -p 12

strawberry ../alignments/align/snowfly_ab_m_Aligned.sortedByCoord.out.bam -o snowfly_ab_m_strawberry.gtf -p 12

strawberry ../alignments/align/snowfly_head_Aligned.sortedByCoord.out.bam -o snowfly_head_strawberry.gtf -p 12

strawberry ../alignments/align/snowfly_ant_Aligned.sortedByCoord.out.bam -o snowfly_ant_strawberry.gtf -p 12

strawberry ../alignments/align/snowfly_thorax_Aligned.sortedByCoord.out.bam -o snowfly_thorax_strawberry.gtf -p 12

# Scallop commands, single-thread
scallop -i ../alignments/align/snowfly_ab_f_Aligned.sortedByCoord.out.bam -o snowfly_ab_f_scallop.gtf

scallop -i ../alignments/align/snowfly_ab_m_Aligned.sortedByCoord.out.bam -o snowfly_ab_m_scallop.gtf

scallop -i ../alignments/align/snowfly_head_Aligned.sortedByCoord.out.bam -o snowfly_head_scallop.gtf

scallop -i ../alignments/align/snowfly_ant_Aligned.sortedByCoord.out.bam -o snowfly_ant_scallop.gtf

scallop -i ../alignments/align/snowfly_thorax_Aligned.sortedByCoord.out.bam -o snowfly_thorax_scallop.gtf