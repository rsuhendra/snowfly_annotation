# Commands used to run Blobtools 

# PREREQUISITES:
# Need to run Blast for taxonomy. Here is the command I used
blastn -db /projects/b1042/RichardSuhendra/nt/nt \
       -query /projects/b1020/Richard/snowfly/data/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads $SLURM_NTASKS_PER_NODE \
       -out blast.out

# You also need to combine all your RNA seq data into one file
# mergebam.txt is just a file with the names of all the files you want to combine
# You need to have a .csi index not a .bai index for blobtools add --cov
# hence the -bc option
samtools merge -b merge_bam.txt --threads 52 snowfly.bam
samtools sort snowfly.bam -o snowfly_sorted.bam -@ 52
samtools index -bc snowfly_sorted.bam
rm snowfly.bam

# Here we use the original, uncleaned assembly fasta file
blobtools create --fasta /projects/b1020/Richard/snowfly/data/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.fasta snowfly

# We add blast hits for taxonomy based on blasting the assembly fasta file
blobtools add --hits blast/blast.out --taxrule bestsumorder --taxdump taxdump snowfly

# Adds coverage data for a specific plot
blobtools add --cov /projects/b1020/Richard/snowfly/align/snowfly_sorted.bam snowfly

# This will run the commands on the cluster, so you can access remotely.
blobtools view --remote snowfly

# "snowfly_blobtools.csv" is one of the taxonomy information
# file you can get from blobtools. Here we remove anything guaranteed 
# to not be Arthropoda
awk -F ',' '{if ($5=="Arthropoda" || $5=="no-hit") print $10}' snowfly_blobtools.csv > kept_scaffolds.txt

# Here we "clean" the assembly
seqtk subseq /projects/b1020/Richard/snowfly/repeats/05_full_out/Chionea.hifiasm.batch2.15x.v210311.-l3-s0.5.asm.p_ctg_HiC.fasta keep2.txt | fold -b80 > Chionea_cleaned_v2.fasta




