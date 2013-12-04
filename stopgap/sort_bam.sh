samtools view -Sb /home/jason/Dropbox/stopgap-measure/$1.sam > /home/jason/Dropbox/stopgap-measure/$1.bam
samtools sort /home/jason/Dropbox/stopgap-measure/$1.bam /home/jason/Dropbox/stopgap-measure/$1_sorted
samtools index /home/jason/Dropbox/stopgap-measure/$1_sorted.bam
rm $1.bam
