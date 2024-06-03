#Run trimmomatic on all paired end files

#!/bin/bash

#set the working directory
cd /mnt/store/sally_rnaseq/data
pwd

#run a for loop to trim all fastq files in data directory
for infile in *_R1_001.fastq.out.gz
do
	base=$(basename ${infile} _R1_001.fastq.out.gz)
	java -jar /mnt/store/sally_rnaseq/jo_QC/jo_trimming/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${infile} ${base}_R2_001.fastq.out.gz /mnt/store/sally_rnaseq/jo_QC/jo_trimming/${base}_R1_001_paired.fastq.out.gz /mnt/store/sally_rnaseq/jo_QC/jo_trimming/${base}_R1_001_unpaired.fastq.out.gz /mnt/store/sally_rnaseq/jo_QC/jo_trimming/${base}_R2_001_paired.fastq.out.gz /mnt/store/sally_rnaseq/jo_QC/jo_trimming/${base}_R2_001_unpaired.fastq.out.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
done

