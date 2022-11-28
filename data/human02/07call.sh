module load freebayes
freebayes -f humanMitoRef.fa --haplotype-length -1 --min-alternate-fraction 0.01 --min-alternate-count 2 --pooled-continuous -p 1 -X -u -i  *s.bam> humanmitoVarsBwaFreeb12dupsMarked.vcf

