


module load anaconda3/personal
source activate Omega
module load macs


for i in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
  for files in $i/*UNIQUE.bam; do
    bname=`basename $files`
    macs2 callpeak -t $files --outdir $i/MACS2output -n ${bname}
  done
done




#Generates a PDF graph for every single .r model in all folders

  for i in $EPHEMERAL/FEBREFRESH/chipseq/chipseq/GSE*; do
    for folders in $i/MACS2FINALOUTPUT; do
      for rmodels in $folders/*UNIQUE.bam_model.r; do
        bname=`basename $rmodels`
        Rscript $folders/${bname}
      done
    done
  done




for files in ./*UNIQUE.bam; do
  bname=`basename $files`
  echo "The file is $bname"
  macs2 callpeak -t $files -c ./SRR099014.fastq.trimmedAligned.sortedByCoord.out.bam_marked_duplicates.bam --outdir MACS2FINALOUTPUT -n ${bname}
done
