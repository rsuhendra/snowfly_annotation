# COMMANDS TO RUN MIKADO
# Most of the commands are from here and are probably better explained: 
# https://mikado.readthedocs.io/en/stable/Tutorial/

# reccomended prerequisites for Mikado

# Portcullis: it's very easy to run. Just do the full option. Afterwards, use filtered
# pass junctions bed file. Here is the command I used to run Portcullis
portcullis full -t 20 /projects/b1020/Richard/snowfly/data/Chionea_cleaned_v2.softmasked.fasta /projects/b1020/Richard/snowfly/alignments/align/snowfly_sorted.bam

# A Swissprot proteins fasta file. Go on NCBI, pick some taxonomic category, 
# e.g plants, diptera, specific species, etc. Then click proteins, and select
# UniprotKB/Swissprot on the left. You can download all of them by not
# selecting any, then "Send to->File->Format=Fasta"


# Scoring files: Look at mikado configure --help for a list of scoring 
# files (under the --scoring section)


# RUNNING MIKADO

# Mikado configure

# If you have an issue where Mikado configure gives you an empty configuration 
# file, you need to downgrade marshmallow and marshmallowdb (look at Mikado
# github issues for more detail)

GENOME=/projects/b1020/Richard/snowfly/data/Chionea_cleaned_v2.softmasked.fasta

# Mikado configure
mikado configure --list config_diptera.txt --reference $GENOME --mode permissive --scoring HISTORIC/insects.yaml --copy-scoring insects.yaml --junctions /projects/b1020/Richard/snowfly/mikado/portcullis/portcullis_out/3-filt/portcullis_filtered.pass.junctions.bed -bt /projects/b1020/Richard/other_data/proteins/uniprot_sprot_diptera/uniprot_sprot_diptera.fasta configuration.yaml


# Mikado prepare
mikado prepare --json-conf configuration.yaml


# Prepare Blast database
makeblastdb -in uniprot_sprot_diptera.fasta -dbtype prot -parse_seqids > blast_prepare.log

# Now we need to run blastx. You NEED to run this command on a cluster or 
# it will be very slow. 
blastx -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" -num_threads $SLURM_NTASKS_PER_NODE -query mikado_prepared.fasta -db /projects/b1020/Richard/other_data/proteins/uniprot_sprot_diptera/uniprot_sprot_diptera.fasta -out mikado_prepared.blast.tsv


# Transdecoder
TransDecoder.LongOrfs -t mikado_prepared.fasta

# This step should also be run on a cluster
blastp -query mikado_prepared.fasta.transdecoder_dir/longest_orfs.pep -db /projects/b1020/Richard/other_data/proteins/uniprot_sprot_diptera/uniprot_sprot_diptera.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads $SLURM_NTASKS_PER_NODE > transdecoder_blastp.outfmt6

TransDecoder.Predict -t mikado_prepared.fasta --retain_blastp_hits transdecoder_blastp.outfmt6

# Mikado serialise
mikado serialise --json-conf configuration.yaml --xml mikado_prepared.blast.tsv --blast_targets /projects/b1020/Richard/other_data/proteins/uniprot_sprot_diptera/uniprot_sprot_diptera.fasta --orfs mikado_prepared.fasta.transdecoder.gff3 --junctions /projects/b1020/Richard/snowfly/mikado/portcullis/portcullis_out/3-filt/portcullis_filtered.pass.junctions.bed


# Mikado pick
mikado pick --configuration configuration.yaml



