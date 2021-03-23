modula load bwa
module load samtools
bwa mem -t 10 -R '@RG\tID:SRR6214421\tSM:SRR6214421' SRR6214434_mt.fa SRR6214421_1.fastq-trimmed-pair1.fastq SRR6214421_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214421j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214423\tSM:SRR6214423' SRR6214434_mt.fa SRR6214423_1.fastq-trimmed-pair1.fastq SRR6214423_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214423j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214425\tSM:SRR6214425' SRR6214434_mt.fa SRR6214425_1.fastq-trimmed-pair1.fastq SRR6214425_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214425j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214426\tSM:SRR6214426' SRR6214434_mt.fa SRR6214426_1.fastq-trimmed-pair1.fastq SRR6214426_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214426j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214429\tSM:SRR6214429' SRR6214434_mt.fa SRR6214429_1.fastq-trimmed-pair1.fastq SRR6214429_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214429j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214431\tSM:SRR6214431' SRR6214434_mt.fa SRR6214431_1.fastq-trimmed-pair1.fastq SRR6214431_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214431j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214433\tSM:SRR6214433' SRR6214434_mt.fa SRR6214433_1.fastq-trimmed-pair1.fastq SRR6214433_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214433j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214434\tSM:SRR6214434' SRR6214434_mt.fa SRR6214434_1.fastq-trimmed-pair1.fastq SRR6214434_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214434j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214436\tSM:SRR6214436' SRR6214434_mt.fa SRR6214436_1.fastq-trimmed-pair1.fastq SRR6214436_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214436j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214438\tSM:SRR6214438' SRR6214434_mt.fa SRR6214438_1.fastq-trimmed-pair1.fastq SRR6214438_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214438j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214440\tSM:SRR6214440' SRR6214434_mt.fa SRR6214440_1.fastq-trimmed-pair1.fastq SRR6214440_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214440j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214420\tSM:SRR6214420' SRR6214434_mt.fa SRR6214420_1.fastq-trimmed-pair1.fastq SRR6214420_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214420j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214422\tSM:SRR6214422' SRR6214434_mt.fa SRR6214422_1.fastq-trimmed-pair1.fastq SRR6214422_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214422j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214424\tSM:SRR6214424' SRR6214434_mt.fa SRR6214424_1.fastq-trimmed-pair1.fastq SRR6214424_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214424j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214427\tSM:SRR6214427' SRR6214434_mt.fa SRR6214427_1.fastq-trimmed-pair1.fastq SRR6214427_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214427j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214428\tSM:SRR6214428' SRR6214434_mt.fa SRR6214428_1.fastq-trimmed-pair1.fastq SRR6214428_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214428j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214430\tSM:SRR6214430' SRR6214434_mt.fa SRR6214430_1.fastq-trimmed-pair1.fastq SRR6214430_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214430j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214432\tSM:SRR6214432' SRR6214434_mt.fa SRR6214432_1.fastq-trimmed-pair1.fastq SRR6214432_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214432j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214435\tSM:SRR6214435' SRR6214434_mt.fa SRR6214435_1.fastq-trimmed-pair1.fastq SRR6214435_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214435j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214437\tSM:SRR6214437' SRR6214434_mt.fa SRR6214437_1.fastq-trimmed-pair1.fastq SRR6214437_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214437j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214439\tSM:SRR6214439' SRR6214434_mt.fa SRR6214439_1.fastq-trimmed-pair1.fastq SRR6214439_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214439j.bam
bwa mem -t 10 -R '@RG\tID:SRR6214441\tSM:SRR6214441' SRR6214434_mt.fa SRR6214441_1.fastq-trimmed-pair1.fastq SRR6214441_1.fastq-trimmed-pair2.fastq | samtools view -b -o SRR6214441j.bam
