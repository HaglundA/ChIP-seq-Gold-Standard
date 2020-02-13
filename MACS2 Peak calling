

module load macs

module load anaconda3/personal
source activate Omega
module load macs


for i in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*; do
  macs2 callpeak -t $i/*duplicates.bam --outdir $i/MACS2output
done
