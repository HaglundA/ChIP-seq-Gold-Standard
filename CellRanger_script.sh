#STEP 1: MAKE GTF FILE

#mkdir refdata-cellranger-mm10-3.0.0_premrna
#cd refdata-cellranger-mm10-3.0.0_premrna

export PATH=/rds/general/user/ah3918/home/cellranger/cellranger-3.1.0:$PATH

wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz

#------Hum

wget ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
gunzip *.gz

#----------------------------------------------------------------------------------------------------------------------


#PBS -lwalltime=23:00:00
#PBS -lselect=1:ncpus=8:ompthreads=8:mem=80gb
export PATH=/rds/general/user/ah3918/home/cellranger/cellranger-3.1.0:$PATH
genomeloc=/rds/general/user/ah3918/ephemeral/cellrangerdata/referenceGenomes/mm10/
cellranger mkgtf $genomeloc/Mus_musculus.GRCm38.93.gtf $genomeloc/Mus_musculus.GRCm38.93.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene

cellranger mkref --genome=mm10 \
                --fasta=$genomeloc/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--genes=$genomeloc/Mus_musculus.GRCm38.93.filtered.gtf \
--ref-version=3.0.0




#PBS -lwalltime=23:00:00
#PBS -lselect=1:ncpus=8:ompthreads=8:mem=80gb

humangenomeloc=/rds/general/user/ah3918/ephemeral/cellrangerdata/referenceGenomes/hg38
export PATH=/rds/general/user/ah3918/home/cellranger/cellranger-3.1.0:$PATH
cellranger mkgtf $humangenomeloc/Homo_sapiens.GRCh38.93.gtf $humangenomeloc/Homo_sapiens.GRCh38.93.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene

cellranger mkref --genome=GRCh38 \
                --fasta=$humangenomeloc/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--genes=$humangenomeloc/Homo_sapiens.GRCh38.93.filtered.gtf \
--ref-version=3.0.0

#----------------------------------------------------------------------------------------------------------------------

#STEP2: CREATE PBSPRO TEMPLATE FILE

#cd cellranger-3.0.2/martian-cs/v3.2.0/jobmanagers/

#cp pbspro.template.example pbspro.template

#PBS -N __MRO_JOB_NAME__
#PBS -V
#PBS -j oe
#PBS -M rf1116@ic.ac.uk
#PBS -l select=1:ncpus=4:mem=60gb
#PBS -l walltime=24:00:00
#PBS -o __MRO_STDOUT__
#PBS -e __MRO_STDERR__

#cd __MRO_JOB_WORKDIR__

#__MRO_CMD__

#STEP3: RUN Cellranger


export PATH=/rds/general/user/ma14618/home/SCRATCH/ax4/CellRanger/cellranger-3.1.0:$PATH

ID=PC
NAME=PC
FQ=/rds/general/user/ma14618/home/SCRATCH/ax4/sourceBiosciencePilot/bp100/HHKHTDRXX/PC

cellranger count --id=$ID --fastqs=$FQ --sample=$NAME --transcriptome=/rds/general/user/ma14618/home/SCRATCH/ax4/CellRanger/refdata-cellranger-mm10-3.0.0_premrna/mm10-3.0.0_premrna.gtf --expect-cells=10000 --jobmode=pbspro --maxjobs=5 --localcores=10 --mempercore=12
