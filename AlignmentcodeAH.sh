#-----------------------------------------------------------------------------
#1. Genome generation and indexing
#-----------------------------------------------------------------------------


module load star/2.7.1a
#MOUSE genome
Genomelocation=/rds/general/user/ah3918/ephemeral/alignments/genomes/mm10
mkdir $Genomelocation/GenomeIndexSTAR


vi STAR.qsub

#Vi text editor
#-----------------------------------------------------------------------------

#PBS -lwalltime=23:00:00
#PBS -lselect=1:ncpus=8:ompthreads=8:mem=80gb
Genomelocation=/rds/general/user/ah3918/ephemeral/alignments/genomes/mm10
module load star/2.7.1a
STAR --runMode genomeGenerate \
--genomeDir $Genomelocation/GenomeIndexSTAR \
--genomeFastaFiles $Genomelocation/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile $Genomelocation/Mus_musculus.GRCm38.98.chr.gtf \
--runThreadN 8

#-----------------------------------------------------------------------------
qsub STAR.qsub

#-----------------------------------------------------------------------------
#2. Alignment
#-----------------------------------------------------------------------------

for folder in ./*; do
  mkdir $folder/alignments
done


GENOME=/rds/general/user/ah3918/ephemeral/alignments/genomes/mm10/GenomeIndexSTAR
ALIGNMENT=/rds/general/user/ah3918/ephemeral/alignments
cd $FILES
ls *.fastq.trimmedx > fastqlist.txt
cd $ALIGNMENT
LIST=/rds/general/user/ah3918/ephemeral/chipseq/GSE69859/fastqlist.txt
name=$(cat $LIST)


#Change GSE directory in code. All GSE folders should have a text file with a list of
#fastq files in the directory


vi STAR.qsub
#-----------------------------------------------------------------------------
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=8:ompthreads=8:mem=96gb
GENOME=/rds/general/user/ah3918/ephemeral/genomes/genomes/mm10/GenomeIndexSTAR
ALIGNMENTSdir=/rds/general/user/ah3918/ephemeral/alignments/FirstAlignment
FILES=/rds/general/user/ah3918/ephemeral/chipseq/GSE106220/*.fastq.trimmed
module load star/2.7.1a

for i in $FILES
do


STAR --genomeDir $GENOME --readFilesIn $i --outFileNamePrefix $ALIGNMENTSdir/ --quantMode GeneCounts --runThreadN 8


done

#-----------------------------------------------------------------------------
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=8:ompthreads=8:mem=96gb
GENOME=/rds/general/user/ah3918/ephemeral/genomes/genomes/mm10/GenomeIndexSTAR
module load star/2.7.1a
for folders in /rds/general/user/ah3918/ephemeral/chipseq/GSE*;do
  for i in $folders/*.fastq.trimmed; do
    STAR --genomeDir $GENOME --readFilesIn $i --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $i  --quantMode GeneCounts --runThreadN 8
  done
done



#MULTIQC
for folders in /rds/general/user/ah3918/ephemeral/Janrefresh/chipseq/GSE*;do
  multiqc $folders
done
