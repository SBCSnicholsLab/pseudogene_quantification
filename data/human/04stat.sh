module load samtools
module load parallel
parallel -j 20 samtools stats {} ">" {.}.stats ::: *d.bam
