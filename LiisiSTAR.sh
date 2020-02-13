# MAPPING WITH STAR
#!/bin/sh 
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=12:mem=60gb
#PBS -M ll3515@ic.ac.uk
#PBS -m bea
#PBS -j oe

# ------------------------------------------------------------------------------------
# ALIGNEMENT and QUANTIFICATION WITH STAR
# ------------------------------------------------------------------------------------
module load star/2.5.0a

# -------------------------------------
# 1. Generate genome (can be skipped if compiled genome already exists)
# -------------------------------------

# I have already done that, but feel free to repeat, if you like!

# # genome fasta and .gtf files for the genome version of interest:
# wherethefilesat=/scratch2/ll3515/Reference_genome_hg38/
# 
# mkdir $wherethefilesat/GenomeIndexSTAR
# 
# STAR --runMode genomeGenerate \
#  --genomeDir $wherethefilesat/GenomeIndexSTAR \
#  --genomeFastaFiles $wherethefilesat/UCSC_hg38/hg38.fa \
#  --sjdbGTFfile $wherethefilesat/UCSC_hg38/hg38.gtf \
#  --runThreadN 12
# 

# -------------------------------------
#2. Alignement
# -------------------------------------

# Specify genome folder (generated in previous step)
GENOME=/scratch2/ll3515/Reference_genome_hg38/GenomeIndexSTAR

FILES="/scratch2/ll3515/DATA/GSE57872/sc_GBM/sample_names.txt"
name=$(cat $FILES)

for i in $name
do

cd /scratch2/ll3515/DATA/GSE57872/sc_GBM/$i


read1=*_1.fastq.gz # fastq file for read1
read2=*_2.fastq.gz # fastq file for read2


STAR --genomeDir $GENOME --readFilesIn $read1 $read2 --readFilesCommand zcat \
 --quantMode GeneCounts --runThreadN 12


done


# ------------------------------------------------------------------------------------
# Run MultiQC
# ------------------------------------------------------------------------------------
cd /scratch2/ll3515/DATA/GSE57872/sc_GBM

multiqc .


# ------------------------------------------------------------------------------------
# Convert SAM to BAM and sort
# ------------------------------------------------------------------------------------
module load samtools/1.3.1

FILES="/scratch2/ll3515/DATA/GSE57872/sc_GBM/sample_names.txt"
name=$(cat $FILES)

for i in $name
do

WRKDIR=/scratch2/ll3515/DATA/GSE57872/sc_GBM/

# Convert SAM to BAM
samtools view -b -S $WRKDIR/$i/Aligned.out.sam > $WRKDIR/$i/Aligned.out.bam

# Sort
samtools sort $WRKDIR/$i/Aligned.out.bam -o $WRKDIR/$i/Aligned.out.sorted.bam

done

# ------------------------------------------------------------------------------------
# DONE! 
# ------------------------------------------------------------------------------------

















