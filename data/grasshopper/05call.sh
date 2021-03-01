#module load freebaes/1.2.0
# freebayes -f ref_11.fasta --haplotype-length -1 --min-alternate-fraction 0.01 --min-alternate-count 2 --pooled-continuous -p 1 -X -u -i  *s.bam> mitoVarsFreeb12.vcf

freebayes -f ref_11.fasta --haplotype-length -1 --min-alternate-fraction 0.01 --min-alternate-count 2 --pooled-continuous -p 1 -X -u -i  *d.bam> mitoVarsFreeb12dupsMarked.vcf

