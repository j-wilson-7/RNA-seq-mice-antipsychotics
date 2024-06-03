#Run FastQC on all the trimmomatic output files

#!/bin/bash
cd /mnt/store/sally_rnaseq/jo_QC/FastQC
pwd

for file in /mnt/store/sally_rnaseq/jo_QC/jo_trimming/data_trimmed_truseq2/paired/*
do
	./fastqc -f fastq ${file}
done