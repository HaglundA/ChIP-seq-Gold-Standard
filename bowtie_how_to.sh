1. ALIGNEMENT
###########################################################################################

#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=30gb
#PBS -M ll3515@ic.ac.uk
#PBS -m bea
#PBS -j oe

# BOWTIE

#module load bowtie/2.1.0 # complains about bowtie2, even though I called just bowtie
module load bowtie/1.0.0

INDEX=/scratch2/ll3515/Reference_genome_hg38/Bowtie_index

cd /scratch2/ll3515/DATA/GSE49402
# <ebwt> base name: GCA_000001405.15_GRCh38_no_alt_analysis_set
bowtie $INDEX/GCA_000001405.15_GRCh38_no_alt_analysis_set -q SRR952608_1_IgG_Rabbit.fastq -m 1 -3 3 -S 2>SRR952608.out > SRR952608.sam

###########################################################################################
# RESULTS: 
SRR952401 (E2F2)
# reads processed: 31817127
# reads with at least one reported alignment: 21919234 (68.89%) <- reads mapped!!!
# reads that failed to align: 461875 (1.45%)
# reads with alignments suppressed due to -m: 9436018 (29.66%)
#Reported 21919234 alignments to 1 output stream(s)

SRR952608 (IgG)
# reads processed: 34143192
# reads with at least one reported alignment: 22532159 (65.99%)
# reads that failed to align: 1164333 (3.41%)
# reads with alignments suppressed due to -m: 10446700 (30.60%)
#Reported 22532159 alignments to 1 output stream(s)


###########################################################################################
2. SAM TO BAM
module load samtools
samtools view -bS SRR952608.sam > SRR952608.bam
###########################################################################################


###########################################################################################
# FAILS TO COMPILE!


2.1 QUALITY CONTROL with PhantomPeakQualTools (ENCODE)
	
	#a) convert the SAM file into a sorted BAM file
#module load samtools
#samtools view -bS SRR952401.sam > SRR952401.bam
	
	# convert the BAM file into TagAlign format, specific to the program that calculates the quality metrics
samtools view -F 0x0204 -o - SRR952401.bam | awk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' | gzip -c > SRR952401.tagAlign.gz
samtools view -F 0x0204 -o - SRR952608.bam | awk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' | gzip -c > SRR952608.tagAlign.gz	
	# run phantompeaktools
	
#Kirill >>>>	
	
PEAKTOOLS=/Data/tools/phantompeaktools/phantompeakqualtools


Rscript ~/phantompeakqualtools/run_spp.R -c=SRR952401.tagAlign.gz -savp -out=SRR952401
 
#Error in library(spp) : there is no package called ‘spp’

###########################################################################################

#>>> Locally >>>
3. CALL PEAKS
# Set up envionment
cd /Users/ll3515/Dropbox/TOOLBOX/MACS2-2.1.1.20160309
export PYTHONPATH=$HOME/lib/python2.7/site-packages/:$PYTHONPATH
export PATH=$HOME/bin/:$PATH
python setup.py install --prefix=$HOME

cd /Users/ll3515/Data/GSE49402/bams
macs2 callpeak --treatment SRR952401.bam --control SRR952608.bam --format BAM -g hs --name E2F2 -B -q 0.01 

# Run the R-script (get plots)
Rscript E2F2_model.r


###########################################################################################


4. GREAT

# E2F2_peaks.narrowPeak >
# Let's format the file as a 3 fields BED file and focus on more significant peaks filtering on q-values.
#awk '$9>40'  E2F2_peaks.narrowPeak | cut -f 1-3 | sed 's/^/chr/' >  E2F2_peaks.narrowPeak.bed # if not UCSC 
awk '$9>40'  E2F2_peaks.narrowPeak | cut -f 1-3  >  E2F2_peaks.narrowPeak.bed


# liftover for Great (locally)
# usage: liftOver oldFile map.chain newFile unMapped

#oldFile and newFile are in bed format by default, but can be in GFF 
#The map.chain file has the old genome as the target and the new genome as the query.

TOOLBOX=/Users/ll3515/Dropbox/TOOLBOX/liftOver
$TOOLBOX/liftOver E2F2_peaks.narrowPeak.bed $TOOLBOX/hg38ToHg19.over.chain E2F2_peaks.hg19.bed E2F2_unmapped.hg19.bed

#inp. to GREAT: http://bejerano.stanford.edu/great/public/html/
















