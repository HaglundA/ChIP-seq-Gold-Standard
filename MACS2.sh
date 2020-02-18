

module load macs
module load anaconda3/personal
source activate Omega
module load macs


for i in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
  for files in $i/*UNIQUE.bam; do
    bname=`basename $files`
    macs2 callpeak -t $files --outdir $i/MACS2output -n ${bname}
  done
done

for i in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
  for folders in $i/MACS2output; do
      Rscript $folders/*UNIQUE.bam_model.r
    done
  done



#Generates a PDF graph for every single .r model in all folders

  for i in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
    for folders in $i/MACS2output; do
      for rmodels in $folders/*UNIQUE.bam_model.r; do
        bname=`basename $rmodels`
        Rscript ${bname}
      done
    done
  done
