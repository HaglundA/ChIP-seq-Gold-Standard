



#= export PATH=/rds/general/user/ah3918/home/JULIA/julia-1.0.5/bin:$PATH


This adds a usermade depositpath (which allows it to be writeable)

export PATH=/rds/general/user/ah3918/home/JULIA/julia-1.0.5/bin:$PATH
export JULIA_DEPOT_PATH=/rds/general/user/ah3918/home/JULIA/julia-1.0.5/PACKAGES:$JULIA_DEPOT_PATH

# " ] " is the Julia command to go into the Pkg REPL


#tutorial
/rds/general/user/ah3918/ephemeral/NETWORKINFERENCE/Tutorial/network_inference_tutorials/simulated_datasets
================ JULIA ===============================================#


#useful julia commands
] test #Tests package installation (in Pkg REPL)
typeof() #R class()
size()   #R dim()
readdir()  #$ ls, R list.files()

; #supresses output



using Pkg
Pkg.add("DelimitedFiles")
Pkg.add("NetworkInference")
Pkg.add("CSV")
Pkg.add("RData")
Pkg.add("LightGraphs")
Pkg.add("Cairo")
Pkg.add("Compose")

using NetworkInference, Cairo, Compose, LightGraphs, GraphPlot, DelimitedFiles

] test NetworkInference

#If OK, proceed

cd("/rds/general/user/ah3918/ephemeral/NETWORKINFERENCE")

using RData
expressionmatrix=load("/rds/general/user/ah3918/ephemeral/SCANALYSIS/HS100bp/Astrocyte_Raw_Expression_Matrix.rds",convert=true)

using CSV
expressionDataFrame=CSV.read("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/NETWORKINFERENCE/Astrocyte_raw_Expression_Matrix.csv")


infer_network(MyDataset,PIDCNetworkInference(),out_file_path="astroDataPreProcessed_Network_17MAR20.txt")
#=============== TUTORIAL ===========================================================#


number_of_genes = 50
organism = "yeast1"
dataset_size = "large"

TutorialDataSet = string("/rds/general/user/ah3918/ephemeral/NETWORKINFERENCE/Tutorial/network_inference_tutorials/simulated_datasets/", number_of_genes, "_", organism, "_", dataset_size, ".txt")

#Choose algorithm

algorithm = PIDCNetworkInference();

#Set threshold. Keep top 85%
threshold = 0.15;

# the @ is to inform Julia simple lang rules might not apply, i.e when calling a function (such as here with 'time'
#This gets the genes and tells you the amount of time it took

@time genes = get_nodes(TutorialDataSet);

network=InferredNetwork(algorithm,genes);



adjacency_matrix, labels_to_ids, ids_to_labels = get_adjacency_matrix(network, threshold);
graph = LightGraphs.SimpleGraphs.SimpleGraph(adjacency_matrix)

#Get Node Labels
number_of_nodes = size(adjacency_matrix)[1]
nodelabels = []
for i in 1 : number_of_nodes
    push!(nodelabels, ids_to_labels[i])
end


#Plot graphs


TutorialPlot=gplot(graph, nodelabel = nodelabels);

using Cairo, Compose

cairo_create

draw(PDF("TutorialTestPlot.png", 10cm, 10cm),TutorialPlot)
