




for files in /rds/general/user/ah3918/ephemeral/FEBREFRESH/chipseq/chipseq/GSE*; do
  for i in $files/*sortedByCoord.out.bam; do
    fastqc -d $files/fastqctemp -o $files/QCAlignedBAM $i
    multiqc $files/QCAlignedBAM -o $files/QCAlignedBAM
  done
  for d in $files/*duplicates.bam; do
    fastqc -d $files/fastqctemp -o $files/QCAlignedNoDupsBAM $d
    multiqc $files/QCAlignedNoDupsBAM -o $files/QCAlignedNoDupsBAM
  done
  for c in $files/*UNIQUE.bam; do
    fastqc -d $files/fastqctemp -o $files/UNIQUEREADSQC $c
    multiqc $files/UNIQUEREADSQC -o $files/UNIQUEREADSQC
  done
done


for files in /rds/general/user/ah3918/ephemeral/FEBREFRESH/chipseq/chipseq/GSE*; do
  mkdir $files FINALQCREPORTS
done



for i in /rdsgpfs/general/ephemeral/user/ah3918/ephemeral/FEBREFRESH/chipseq/chipseq/GSE106220; do
  multiqc $i/UNIQUEREADSQC $i/QCAlignedNoDupsBAM $i/QCAlignedBAM -o $i/FINALQCREPORTS -n "Finalreport"
done
