# Snowfly Annotation

This directory contains scripts used to run the pipeline from *chionea et al*. Here, we succinctly describe the use, and order for each script. A broad overview can be obtained in the methods,and more specific instructions can be found in each shell script. 

### Step 1: Align reads to the genome. 
Refer to *align.sh*

### Step 2: Run Blobtools to clean the assembly if there is some contamination
Refer to *blobtools.sh*

### Step 3: Realign if scaffolds were removed in Step 2
Again refer to *align.sh*

### Step 4: Run RepeatModeler and RepeatMasker
Refer to *repeatmodeler.sh* and *repeatmasker.sh*

### Step 5: Run Braker
Refer to *braker.sh*

### Step 6: Run transcriptome tools Stringtie, Strawberry, and Scallop
Refer to *transcriptomes.sh*

### Step 7: Run Mikado to combine all types and generate the final annotation
Refer to  *mikado.sh* and *config.txt*