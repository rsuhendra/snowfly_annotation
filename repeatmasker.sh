#!/bin/bash
#SBATCH -A XXXXX
#SBATCH -p short
#SBATCH --job-name="snowfly repeatmasker"
#SBATCH -N 1
#SBATCH --ntasks-per-node=52
#SBATCH -t 8:00:00
#SBATCH --mem=100G
#SBATCH --output=repeatmasker.log

cd $SLURM_SUBMIT_DIR

module purge all
source activate augusco

# Here we run RepeatMasker in rounds for better understanding. Firstly, simple repeats 
# are singled out. Then, we use mask repeats found in the repbase library. Then, 
# we mask using our results from RepeatModeler (seperated into known and unknown parts).
# Finally, we combine our the results and softmask the genome as is required for BRAKER.

# RepeatMasker Round 1
RepeatMasker -pa 13 -a -e ncbi -dir 01_simple_out -noint -xsmall /projects/b1020/Richard/snowfly/data/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.fasta

rename fasta simple_mask 01_simple_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC*
rename .masked .masked.fasta 01_simple_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC*

# RepeatMasker Round 2
RepeatMasker -pa 13 -a -e ncbi -dir 02_repbase_out -nolow -species Diptera /projects/b1020/Richard/snowfly/repeats/01_simple_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.simple_mask.masked.fasta

rename simple_mask.masked.fasta repbase_mask 02_repbase_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC*
rename .masked .masked.fasta 02_repbase_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC*

# RepeatMasker Round 3
RepeatMasker -pa 13 -a -e ncbi -dir 03_known_out -nolow -lib /projects/b1020/Richard/snowfly/repeats/snowfly-families.prefix.fa.known /projects/b1020/Richard/snowfly/repeats/02_repbase_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.repbase_mask.masked.fasta

rename repbase_mask.masked.fasta known_mask 03_known_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC*
rename .masked .masked.fasta 03_known_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC*

# RepeatMasker Round 4
RepeatMasker -pa 13 -a -e ncbi -dir 04_unknown_out -nolow -lib /projects/b1020/Richard/snowfly/repeats/snowfly-families.prefix.fa.unknown /projects/b1020/Richard/snowfly/repeats/03_known_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.known_mask.masked.fasta

rename known_mask.masked.fasta unknown_mask 04_unknown_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC*
rename .masked .masked.fasta 04_unknown_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC*


# combine tabular files
cat 01_simple_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.simple_mask.out <(cat 02_repbase_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.repbase_mask.out | tail -n +4) <(cat 03_known_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.known_mask.out | tail -n +4) <(cat 04_unknown_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.unknown_mask.out | tail -n +4) > 05_full_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.full_mask.out


# Convert out files to gff3
/projects/b1020/Richard/snowfly/repeats/rmOutToGFF3custom -o 05_full_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.full_mask.out > 05_full_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.full_mask.gff3


# Soft mask the genome 
bedtools maskfasta -fi /projects/b1020/Richard/snowfly/data/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.fasta -bed /projects/b1020/Richard/snowfly/repeats/05_full_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.full_mask.gff3 -fo Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.softmasked.fasta -soft
