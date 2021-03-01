module load sratools
while read line ; do fasterq-dump -p $line; done < sampleList
