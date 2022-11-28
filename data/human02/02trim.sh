module load skewer
module load pigz
for i in *gz; do skewer -q 30 $i; done

module load fastqc
fastqc -t 10 *fastq
pigz -p 10 *fastq
