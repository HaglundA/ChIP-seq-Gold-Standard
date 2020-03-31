

#Installation

module load anaconda3/personal
source activate Omega #Python
conda install wget samtools r-essentials bioconductor-deseq2 bioconductor-edger

mkdir HOMER
cd HOMER
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install

export PATH=$PATH:/rds/general/user/ah3918/home/HOMER/bin/



#Configuration - this shows a list of genomes, promoters, and organisms.
#Can be customized as needed

perl /rds/general/user/ah3918/home/HOMER/configureHomer.pl -list






#=========== MOTIF ANALYSIS ===========================#

for folders in /rds/general/user/ah3918/ephemeral/MARREFRESH/chipseq/GSE*; do
  mkdir $folders/HOMEROUTPUT
done

for folders in /rds/general/user/ah3918/ephemeral/MARREFRESH/chipseq/GSE*; do
  for files in $folders/BLACKLISTFILTERED; do
    for beds in $files/*.filtered.bed; do
      bname=`basename $beds`
    findMotifsGenome.pl $beds mm10 $folders/HOMEROUTPUT/$bname
    done
  done
done


#Screen session 14th March 26813.HOMER
#just run screen -r HOMER

findMotifsGenome.pl
