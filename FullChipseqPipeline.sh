

# This is the full pipeline I have developed (from past papers/inhouse code/online recommendations)
# and used for the analysis of public transcription factor ChIP-seq (TF-ChIP) analysis,
# from initial download (SRA) all the way to enrichment.
# A diagram of the general pipeline will soon be featured in the readme.
#
# For any questions/suggestions/comments outside of github, please email me at
# a.haglund@outlook.com
#
# ------------------------------------------------------------------------------
# 1. DOWNLOADING DATASETS (GEO)
# ------------------------------------------------------------------------------
# 1. First, identify GSE dataset (ex. GSE23219)
# 2. Identify the SRA database
# 3. Change numbers accordingly (So far I have found no better way of downloading
# multiple SRA files). Download SRA toolkit first.


export SRA_PATH=yourdirectory/SRAtoolkit/bin

$PATH_SRA/fastq-dump cat sralist.text

#Change "i" according to the last digits of the SRA accession codes to download iteratively

PATH_SRA=/rds/general/user/ah3918/home/sratoolkit/bin
for (( i = 488; i <= 501; i++ ))
  do
  $PATH_SRA/fastq-dump SRR2062$i
done


#FastQC analysis

mkdir yourQCdirectory/FastQCtemp
mkdir yourQCdirectory/FastQCoutput
fastqc -d FastQCtemp -o FastQCoutput *.fastq #<-input

module load multiqc
multiqc ./FastQCoutput
# ------------------------------------------------------------------------------
# 2. TRIMMING (Trimmomatic)
# ------------------------------------------------------------------------------
#Trimmomatic will trim out adapters used during sequencing process.

Trimmomaticpath=yourdirectory/Trimmomatic-0.39
adapterpath=yourdirectory/Trimmomatic-0.39/adapters

#single set trimmomatic

cd yourGSEdirectory
for i in ../*
do trimmomatic SE -phred33 $i.fastq $i.fastqILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done


for i in ../*.fastq
do java -jar $Trimmomaticpath/trimmomatic-0.39.jar SE -phred33 $i $i.trimmed ILLUMINACLIP:$adapterpath/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done


#Nested loop for trimming fastq files in multiple directories

Trimmomaticpath=yourdirectory/Trimmomatic-0.39
adapterpath=yourdirectory/Trimmomatic-0.39/adapters

for folder in /yourdatadirectory/GSE*; do
  for i in $folder/*.fastq; do
    java -jar $Trimmomaticpath/trimmomatic-0.39.jar SE \
    -threads 8 \
    -phred33 $i $i.trimmed ILLUMINACLIP:$adapterpath/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
  done
done

#At this point, FastQC + MultiQC in step1 can be repeated to compare the quality of the reads pre/post-trimming

# ------------------------------------------------------------------------------
# 3. ALIGNMENT/MAPPING (STAR)
# ------------------------------------------------------------------------------

# 1. First, download appropriate genome files (fasta and gtf files). I have been using ENSEMBL mouse genome (mm10)
# 2. Download STAR

module load star/2.7.1a
Genomelocation=yourgenomelocation/mm10
mkdir $Genomelocation/GenomeIndexSTAR

# 3. Generate reference genome
# If using an HPC, this is where i submit a PBS script


vi STARgenome.qsub


# ------------------------------------------------------------------------------#VI TEXT editor
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=8:ompthreads=8:mem=96gb

module load star/2.7.1a
STAR --runMode genomeGenerate \
--genomeDir $Genomelocation/GenomeIndexSTAR \
--genomeFastaFiles $Genomelocation/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile $Genomelocation/Mus_musculus.GRCm38.98.chr.gtf \
--runThreadN 8

# ------------------------------------------------------------------------------

vi STARalign.qsub


# ------------------------------------------------------------------------------#VI TEXT editor

# 3. Align reads

#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=8:ompthreads=8:mem=96gb
GENOME=yourgenomelocation/mm10/GenomeIndexSTAR
module load star/2.7.1a
for folders in chipseq/GSE*;do
  for i in $folders/*.fastq.trimmed; do
    STAR --genomeDir $GENOME --readFilesIn $i --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $i  --quantMode GeneCounts --runThreadN 8
  done
done

# ------------------------------------------------------------------------------

The alignments can then be multiQCd

module load multiqc
for folders in chipseq/GSE*;do
  multiqc $folders
done


# ------------------------------------------------------------------------------
# 4. DUPLICATE REMOVAL/UNIQUE READS FILTERING
# ------------------------------------------------------------------------------
