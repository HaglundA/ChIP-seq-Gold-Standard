


export PATH=/rds/general/user/ah3918/home/JULIA/julia-1.0.5/bin:$PATH
export JULIA_DEPOT_PATH=/rds/general/user/ah3918/home/JULIA/julia-1.0.5/PACKAGES:$JULIA_DEPOT_PATH


using NetworkInference, Cairo, Compose, LightGraphs, GraphPlot, DelimitedFiles


MyDataset=string("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/NETWORKINFERENCE/18MARCH20_Pyramidal_Raw_Expression_MatrixDataPreProcessed.txt")

infer_network(MyDataset,PIDCNetworkInference(),out_file_path="23MAR20_Pyramidal_Network_File.txt")


#23MAR20 - screen 17711
