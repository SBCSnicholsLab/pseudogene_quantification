#module load samtools
for i in *j.bam; do samtools sort -@ 20 -o ${i%.*}s.bam $i; done

