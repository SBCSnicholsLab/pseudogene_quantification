#module load parallel
#module load picard
#parallel -j 5 picard MarkDuplicates I={} O={.}d.bam 
M={.}-dup-metrics.txt ::: *s.bam
