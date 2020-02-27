#Conda environment set up

module load anaconda3/personal
#if first set up
anaconda-setup

#template
conda create -n envname python=2.7

#Environments to date
conda create -n Epsilon julia #Julia Environment
conda create -n Omega python #Python Environment
conda create -n Gamma r #R environment



conda install -c r colorout
conda install -c r dropletutils
conda install -c r seurat
conda install -c r ggplot2


conda install -c r DoubletFinder

conda install -c r dplyr

conda install -c r conos
