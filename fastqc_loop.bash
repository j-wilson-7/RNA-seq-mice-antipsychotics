#Run FastQC on all the files in the data folder

#!/bin/bash
cd /mnt/store/sally_rnaseq/jo_fastQC/FastQC
pwd

for file in /mnt/store/sally_rnaseq/data/*
do
	./fastqc -f fastq ${file}
done