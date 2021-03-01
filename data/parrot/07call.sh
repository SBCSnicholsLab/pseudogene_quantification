module load freebayes/1.2.0

freebayes -f SRR6214434_mt.fa --haplotype-length -1 --min-alternate-fraction 0.01 --min-alternate-count 2 --pooled-continuous -p 1 -X -u -i  *d.bam> parrotMitoVarsBwaFreeb12dupsMarked.vcf

