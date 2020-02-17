
This is all of the code necessary for the generation of a gold standard from public ChIP-seq data.


General Pipeline

1. Download GEO datasets using SRAtoolkit
2. Perform QC (FastQC + MultiQC) to visualize quality of reads
3. Trim adapters, perform QC again, compare
4. Align/Map reads (STAR)
5. Remove duplicates (picardtools) and filter uniquely mapped reads (samtools)
6. Call peaks (MACS2)
