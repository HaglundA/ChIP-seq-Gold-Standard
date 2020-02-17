PATH_SRA=/rds/general/user/ah3918/home/sratoolkit/bin
$PATH_SRA/fastq-dump SRR6223102 SRR6223103 SRR6223104 SRR6223105


$PATH_SRA/fastq-dump cat sralist.text


PATH_SRA=/rds/general/user/ah3918/home/sratoolkit/bin
for (( i = 488; i <= 501; i++ ))
  do
  $PATH_SRA/fastq-dump SRR2062$i
done

#FASTQ QC control
#Set wd as GSE folder
#Make output and temp directory for FastQC


#FastQC analysis

mkdir fastqctemp
mkdir FastQCoutputTrimmed
fastqc -d fastqctemp -o FastQCoutput *.fastq.trimmed
