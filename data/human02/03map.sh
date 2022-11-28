module load bwa
module load samtools
bwa mem -t 10 -R '@RG\tID:ERR251013\tSM:ERR251013' humanMitoRef.fa ERR251013_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR251013j.bam
bwa mem -t 10 -R '@RG\tID:ERR015526\tSM:ERR015526' humanMitoRef.fa ERR015526_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR015526j.bam
bwa mem -t 10 -R '@RG\tID:SRR101474\tSM:SRR101474' humanMitoRef.fa SRR101474_1.fastq-trimmed.fastq.gz | samtools view -b -o SRR101474j.bam
bwa mem -t 10 -R '@RG\tID:ERR013101\tSM:ERR013101' humanMitoRef.fa ERR013101_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR013101j.bam
bwa mem -t 10 -R '@RG\tID:ERR015522\tSM:ERR015522' humanMitoRef.fa ERR015522_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR015522j.bam
bwa mem -t 10 -R '@RG\tID:ERR015525\tSM:ERR015525' humanMitoRef.fa ERR015525_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR015525j.bam
bwa mem -t 10 -R '@RG\tID:ERR162841\tSM:ERR162841' humanMitoRef.fa ERR162841_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR162841j.bam
bwa mem -t 10 -R '@RG\tID:ERR013103\tSM:ERR013103' humanMitoRef.fa ERR013103_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR013103j.bam
bwa mem -t 10 -R '@RG\tID:ERR015531\tSM:ERR015531' humanMitoRef.fa ERR015531_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR015531j.bam
bwa mem -t 10 -R '@RG\tID:ERR018423\tSM:ERR018423' humanMitoRef.fa ERR018423_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR018423j.bam
bwa mem -t 10 -R '@RG\tID:SRR043397\tSM:SRR043397' humanMitoRef.fa SRR043397_1.fastq-trimmed.fastq.gz | samtools view -b -o SRR043397j.bam
bwa mem -t 10 -R '@RG\tID:ERR019907\tSM:ERR019907' humanMitoRef.fa ERR019907_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR019907j.bam
bwa mem -t 10 -R '@RG\tID:ERR016162\tSM:ERR016162' humanMitoRef.fa ERR016162_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR016162j.bam
bwa mem -t 10 -R '@RG\tID:ERR015484\tSM:ERR015484' humanMitoRef.fa ERR015484_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR015484j.bam
bwa mem -t 10 -R '@RG\tID:ERR018436\tSM:ERR018436' humanMitoRef.fa ERR018436_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR018436j.bam
bwa mem -t 10 -R '@RG\tID:ERR016118\tSM:ERR016118' humanMitoRef.fa ERR016118_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR016118j.bam
bwa mem -t 10 -R '@RG\tID:ERR016137\tSM:ERR016137' humanMitoRef.fa ERR016137_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR016137j.bam
bwa mem -t 10 -R '@RG\tID:ERR016064\tSM:ERR016064' humanMitoRef.fa ERR016064_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR016064j.bam
bwa mem -t 10 -R '@RG\tID:ERR016081\tSM:ERR016081' humanMitoRef.fa ERR016081_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR016081j.bam
bwa mem -t 10 -R '@RG\tID:ERR018448\tSM:ERR018448' humanMitoRef.fa ERR018448_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR018448j.bam
bwa mem -t 10 -R '@RG\tID:SRR792951\tSM:SRR792951' humanMitoRef.fa SRR792951_1.fastq-trimmed.fastq.gz | samtools view -b -o SRR792951j.bam
bwa mem -t 10 -R '@RG\tID:ERR013120\tSM:ERR013120' humanMitoRef.fa ERR013120_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR013120j.bam
bwa mem -t 10 -R '@RG\tID:ERR016251\tSM:ERR016251' humanMitoRef.fa ERR016251_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR016251j.bam
bwa mem -t 10 -R '@RG\tID:ERR018469\tSM:ERR018469' humanMitoRef.fa ERR018469_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR018469j.bam
bwa mem -t 10 -R '@RG\tID:ERR016258\tSM:ERR016258' humanMitoRef.fa ERR016258_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR016258j.bam
bwa mem -t 10 -R '@RG\tID:ERR016269\tSM:ERR016269' humanMitoRef.fa ERR016269_1.fastq-trimmed.fastq.gz | samtools view -b -o ERR016269j.bam

