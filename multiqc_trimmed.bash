#Run multiQC on all FastQC files for trimmed data

#!/bin/bash
cd /mnt/store/sally_rnaseq/jo_QC/jo_multiQC/MultiQC
pwd

multiqc /mnt/store/sally_rnaseq/jo_QC/jo_trimming/fastQC_trimmed/paired -o /mnt/store/sally_rnaseq/jo_QC/jo_trimming/multiQC_trimmed
