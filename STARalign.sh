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

#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=8:ompthreads=8:mem=96gb
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
