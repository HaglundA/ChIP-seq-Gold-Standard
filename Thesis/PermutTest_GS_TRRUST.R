options(stringsAsFactors=F)
library(colorout)


# Load networks
load("~/PIDC_networks.RData")
#[1,] "PIDC.NW.list"
str(PIDC.NW.list)

# Function to run permutation
perm.genes <-function(x){
  smp = sample(AllGenes, size = length(TF), replace = F)
  a = NW[NW$V1 %in% smp & NW$V2 %in% smp,]
  out = sum(a$V3)
  return(out)
}


### --------- Permutation GS --------------------------------------------------------------------------

# Load Gold Standard
GS = readRDS("~/GoldStandard.rds")
str(GS)

# for TF that have more than 10K targets, take top 25%

ResGS = list()

for (i in 1:length(PIDC.NW.list)){
  i = 5
  obsList = list()
  diffList = list()
  pvalList = list()
  print(i)

  for (j in 1:length(GS)){
    # CellType
    print(paste("TF nr ",j))
    CellType = names(PIDC.NW.list)[i]

    # Network file
    NW = PIDC.NW.list[[i]]
    AllGenes = unique(c(NW$V1, NW$V2)); #Head(AllGenes)

    # TF targets
    nameTF = names(GS[j])
    TF = GS[[j]]

    # TF targets expressed in cell-type
    TF = TF[TF%in%AllGenes]; # Head(TF)

    # PIDC scores between targets
    mi_TF = NW[NW$V1 %in% TF & NW$V2 %in% TF,]; # Head(mi_TF)

    # The observed value in the test set
    obs=sum(mi_TF$V3); obsList[[nameTF]]<-obs

    # Sequence of randomised results
    diff <- replicate(1000, perm.genes(length(TF))); diffList[[nameTF]]<-diff
    #diff <- "test"; diffList[[nameTF]]<-diff

    # P-value
    res=(sum(abs(diff) > abs(obs)) + 1) / (length(diff) + 1); pvalList[[nameTF]]<-res
    #res="test"; pvalList[[nameTF]]<-res

    ResGS[[CellType]] = list(obs = obsList, diff = diffList, pval = pvalList)

  }
}

str(ResGS)

#save(ResGS, file = "GS_permut_enrich.RData")



### --------- Permutation TRRUST --------------------------------------------------------------------------

# Load publicly available TRRUST TF-target sets
Load("Alexander_Haglund/Project/GOLD_STANDARD/Public_targets/TRRUST_mm.RData")
str(TRRUST)

ResTRRUST = list()

for (i in 1:length(PIDC.NW.list)){
  obsList = list()
  diffList = list()
  pvalList = list()
  print(i)

  for (j in 1:length(TRRUST)){
    # CellType
    print(paste("TF nr ",j))
    CellType = names(PIDC.NW.list)[i]

    # Network file
    NW = PIDC.NW.list[[i]]
    AllGenes = unique(c(NW$V1, NW$V2)); #Head(AllGenes)

    # TF targets
    nameTF = names(TRRUST[j])
    TF = TRRUST[[j]]

    # TF targets expressed in cell-type
    TF = TF[TF%in%AllGenes]; # Head(TF)

    # PIDC scores between targets
    mi_TF = NW[NW$V1 %in% TF & NW$V2 %in% TF,]; # Head(mi_TF)

    # The observed value in the test set
    obs=sum(mi_TF$V3); obsList[[nameTF]]<-obs

    # Sequence of randomised results
    diff <- replicate(1000, perm.genes(length(TF))); diffList[[nameTF]]<-diff
    #diff <- "test"; diffList[[nameTF]]<-diff

    # P-value
    res=(sum(abs(diff) > abs(obs)) + 1) / (length(diff) + 1); pvalList[[nameTF]]<-res
    #res="test"; pvalList[[nameTF]]<-res

    ResTRRUST[[CellType]] = list(obs = obsList, diff = diffList, pval = pvalList)

  }
}

str(ResTRRUST)

#save(ResTRRUST, file = "TRRUST_permut_enrich.RData")
