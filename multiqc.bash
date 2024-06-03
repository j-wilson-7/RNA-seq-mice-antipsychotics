#Run multiQC on all FastQC files in a certain folder

#!/bin/bash
cd /mnt/store/sally_rnaseq/jo_multiQC/MultiQC
pwd

multiqc /mnt/store/sally_rnaseq/data -o /mnt/store/sally_rnaseq/jo_multiQC
