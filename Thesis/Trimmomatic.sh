Trimmomaticpath=/rds/general/user/ah3918/home/Trimmomatic-0.39
adapterpath=/rds/general/user/ah3918/home/Trimmomatic-0.39/adapters

#single set trimmomatic
cd yourGSEdirectory
for i in ../*
do trimmomatic SE -phred33 $i.fastq $i.fastqILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done


for i in ../*.fastq
do java -jar $Trimmomaticpath/trimmomatic-0.39.jar SE -phred33 $i $i.trimmed ILLUMINACLIP:$adapterpath/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done
#-------------------------------------------
#Nested loop for trimming fastq files in multiple directories (assuming you are in the directory where all those datasets are


#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=8:ompthreads=8:mem=96gb
Trimmomaticpath=/rds/general/user/ah3918/home/Trimmomatic-0.39
adapterpath=/rds/general/user/ah3918/home/Trimmomatic-0.39/adapters
for folder in /rds/general/user/ah3918/ephemeral/chipseq/*; do
  for i in $folder/*.fastq; do
    java -jar $Trimmomaticpath/trimmomatic-0.39.jar SE \
    -threads 8 \
    -phred33 $i $i.trimmed ILLUMINACLIP:$adapterpath/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
  done
done

#---------------------------------
