source("~/000.R_functions.R")

###  Analyse inferred PIDC networks

# - use chosen threshold to create igraph graph for filtered and imputed data
# - cluster to identify gene modules
# - compare cluster similarities for different datasets / clustering algorithms
# - GO analysis of larger clusters


### --------- Load packages ----------------------------------------------------------------------
# NB note that some of these packages mask functions from other packages (e.g. topGO masks some of igraph)
# so use namespace to refer to functions if needed, e.g. igraph::sizes(x)

library(igraph)             # for graph analysis
library(fields)             # for easy heatmap plot with legend (image.plot)
library(data.table)         # for fast file import (fread)
library(biomaRt)            # for ensemblID - gene name conversions
library(topGO)              # for GO analysis
library(org.Hs.eg.db)       # for GO annotations
library(gdata)              # for writing GO table results to file
library(RColorBrewer)       # for colours
library(spaa)               # for matrix to list

source('~/functions/3_analyse_networks_fns.R')      # associated functions
source('~/4_plot_graphs_fns.R')      # associated functions


### --------- Set variables and directories --------------------------------------------------------

# working directory
wd <-"~/NETWORKINFERENCE/"

# variables
threshold_pc = 0.0002        # threshold for network edges (selects top X*100 % of edges)
min_size = 30               # min size of clusters to analyse (in comparisons, GO analysis etc)
cluster_alg = 'ml'          # clustering algorithm (ml = louvain)

# dir
dir_in = paste0(wd,'network_processing/input/')          # location of PIDC output txt files
dir_out = paste0(wd,'network_processing/', 'pc', threshold_pc*100, '/')          # location to save output files
dir_fig = paste0(wd,'network_processing/figures/', 'pc', threshold_pc*100, '/')          # location to save figures

if(!file.exists(dir_out)){dir.create(dir_out)}  # create directories if absent
if(!file.exists(dir_fig)){dir.create(dir_fig)}



### --------- Import data --------------------------------------------------------------------------

# import data (pidc network inference results) from input directory
pidc_results = load_data(dir_in)
networks = pidc_results$pidc_results            # list with entry per inferred network, list of edges in PIDC score order
nodes = pidc_results$nodes                      # reference list of genes in alphabetical order (same for all networks)
rm(pidc_results)
names(networks) = c("oligodendrocytes")

### --------- Create igraphs -----------------------------------------------------------------------

# create igraphs using chosen threshold (and record membership of largest component for each graph)
networks_g = get_thresholded_graph(networks, nodes, threshold_pc, dir_out)
list_graphs = networks_g$list_graphs            # list with entry (igraph) per network, created using chosen threshold
node_comp = networks_g$node_comp                # matrix indicating membership of largest component for each node in nodes
rm(networks_g)



### --------- Plot igraphs -----------------------------------------------------------------------

# threshold to get list of edges from chosen graphs
g_adj = get_adjacency_matrix(networks[[1]], nodes, threshold_pc)

# itrate through different thrsholds to
#p = c(0.0002, 0.0005, 0.001, 0.002)
p = c(0.002, 0.005, 0.01)

meanConnectivity=list()
medianConnectivity=list()

for(iter in p){

threshold_pc = iter
g_adj = get_adjacency_matrix(networks[[1]], nodes, threshold_pc)
rownames(g_adj)<-nodes; colnames(g_adj)<-nodes

pairWise = dist2list(as.dist(t(g_adj)))
pairWise = pairWise[pairWise$value == 1,]

temp = as.data.table(pairWise)
edges = temp[seq(2, dim(temp)[1], 2)]             # remove duplicates
pairWise = as.data.frame(edges)

unique_genes <- unique(c(as.character(pairWise[,"col"]),as.character(pairWise[,"row"]))); Head(unique_genes)

GeneConnectivity = list()
for(i in unique_genes){
  GeneConnectivity[[i]]<-nrow(pairWise[pairWise[,"col"] == i | pairWise[,"row"] == i,])
}
GeneConnectivity = sort(unlist(GeneConnectivity), decreasing =T); head(GeneConnectivity)
#Gpc5  Lsamp   Apoe  Npas3    Ntm Slc1a2
#  30     23     22     14     13     13
mean(GeneConnectivity) # 2.54251
median(GeneConnectivity)

meanConnectivity[[paste0("t_",iter)]]<-mean(GeneConnectivity)
medianConnectivity[[paste0("t_",iter)]]<-median(GeneConnectivity)

}
print(medianConnectivity)
print(meanConnectivity)

# choose a thershold that results in mean connectivity between 2 and 5

threshold_pc = 0.002

g_adj = get_adjacency_matrix(networks[[1]], nodes, threshold_pc)
rownames(g_adj)<-nodes; colnames(g_adj)<-nodes

pairWise = dist2list(as.dist(t(g_adj)))
pairWise = pairWise[pairWise$value == 1,]

temp = as.data.table(pairWise)
edges = temp[seq(2, dim(temp)[1], 2)]             # remove duplicates
pairWise = as.data.frame(edges)
Head(pairWise)



# For CytoScape
write.delim(pairWise, file = "~/Oligodendrocytes_pc0.2.txt", row.names = F)



# get graph layout
#g = graph_from_adjacency_matrix(g_adj, mode = 'undirected')
#edge_colours = adjustcolor('red', alpha = 0.1)
#l = layout_with_fr(g)
#l = layout_with_kk(g)
#l = layout_with_lgl(g)

#pdf("test_plot.pdf")
#plot(g, layout = l,
#     vertex.label = '', vertex.size = 2,
#     vertex.frame.color = adjustcolor('gray50', alpha = 0.4),
#     vertex.color = temp_x,
#     edge.color = edge_colours
#)
#dev.off()
