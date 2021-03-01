#module load samtools
ls *d.bam | parallel -j 20 samtools stats {} ">" {.}.stats ::: *d.bam
