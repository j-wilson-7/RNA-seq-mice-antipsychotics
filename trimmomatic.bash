#Run trimmomatic on specific forward and reverse paired ends

#!/bin/bash

#set the working directory
cd /mnt/store/sally_rnaseq/data
#cd /mnt/store/sally_rnaseq/jo_QC/jo_trimming/Trimmomatic-0.39
pwd

#run trimmomatic-0.39 paired end mode
java -jar /mnt/store/sally_rnaseq/jo_QC/jo_trimming/Trimmomatic-0.39/trimmomatic-0.39.jar PE CLOZ2_GTCCGC_L004_R1_001.fastq.out.gz CLOZ2_GTCCGC_L004_R2_001.fastq.out.gz /mnt/store/sally_rnaseq/jo_QC/jo_trimming/CLOZ2_GTCCGC_L004_R1_001_paired.fastq.out.gz /mnt/store/sally_rnaseq/jo_QC/jo_trimming/CLOZ2_GTCCGC_L004_R1_001_unpaired.fastq.out.gz /mnt/store/sally_rnaseq/jo_QC/jo_trimming/CLOZ2_GTCCGC_L004_R2_001_paired.fastq.out.gz /mnt/store/sally_rnaseq/jo_QC/jo_trimming/CLOZ2_GTCCGC_L004_R2_001_unpaired.fastq.out.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
