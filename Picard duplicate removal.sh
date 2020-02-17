

#Nested loop to remove all duplicates in every #GSE subfolder and store them in files called *_marked_duplicates.bam

for i in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
  for file in $i/*.bam; do
    bname=`basename $file`
      echo "The bam file is:" $bname
      sample=$(basename $file .bam | cut -d- -f1)
      picard MarkDuplicates \
      I=$i/$bname \
      O=$i/${bname}_marked_duplicates.bam \
      M=$i/${bname}_marked_dup_metrics.txt \
      READ_NAME_REGEX=null \
      REMOVE_DUPLICATES=true
    done
done

#QC every single bam file including duplicate removal bam files.
#This is because MACS2 apparently removes duplicates, so I want
#to run MACS2 on both dupremoved and dupkept, but first compare QC



for i in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
  mkdir "$i"/QCAlignedBAM
  mkdir "$i"/QCAlignedNoDupsBAM
done

module load anaconda3/personal
source activate Omega
module load fastqc
module load multiqc
for files in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
  for i in $files/*sortedByCoord.out.bam; do
    fastqc -d $files/fastqctemp -o $files/QCAlignedBAM $i
    multiqc $files/QCAlignedBAM -o $files/QCAlignedBAM
  done
  for d in $files/*duplicates.bam; do
    fastqc -d $files/fastqctemp -o $files/QCAlignedNoDupsBAM $d
    multiqc $files/QCAlignedNoDupsBAM -o $files/QCAlignedNoDupsBAM
  done
done
for files in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
  for d in $files/*UNIQUE.bam; do
    fastqc -d $files/fastqctemp -o $files/UNIQUEREADSQC $d
    multiqc $files/UNIQUEREADSQC -o $files/UNIQUEREADSQC
done



module load samtools
for i in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
  for files in $i/*duplicates.bam; do
    samtools view -q 4 -b $files > $i/UNIQUEREADS/"$files"UNIQUE.bam
  done
done


fastqc -d fastqctemp -o FastQCoutput ../SRR6223105.fastq.trimmedAligned.sortedByCoord.out.bam_marked_duplicates.bam
