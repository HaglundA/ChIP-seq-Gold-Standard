#Gamma2

library(ChIPseeker)
library(GenomicRanges)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

SRRList<-read.delim("/rds/general/user/ah3918/ephemeral/MARREFRESH/GRANGES/SRRList.txt")
setwd("/rds/general/user/ah3918/ephemeral/MARREFRESH/GRANGES/PeakFilesOnly")


fileNames <-list.files()
SampleNames <-gsub("_peaks.xls","",fileNames)

macsPeaks_DF_list <- list()

for(i in 1:length(SampleNames)){
  macsPeaks_DF_list[[SampleNames[[i]]]] <- read.delim(fileNames[i], comment.char="#", header = T, check.names=F)
  print(i)
}


GrangeObjectsList<-list()



for(i in 1:length(macsPeaks_DF_list)){
  GrangeObjectsList[[i]]<-GRanges(
   seqnames=as.data.frame(macsPeaks_DF_list[[i]])[,"chr"],
   IRanges(as.data.frame(macsPeaks_DF_list[[i]])[,"start"],
          (as.data.frame(macsPeaks_DF_list[[i]])[,"end"])
        )
      )
}



# add summit position and fold enrichment
for(i in 1:length(macsPeaks_DF_list)){
mcols(GrangeObjectsList[[i]]) <- as.data.frame(macsPeaks_DF_list[]
  [[i]])[,c("abs_summit", "fold_enrichment")]
}


GrangeObjectsList[[1]][1:2,]

for(i in 1:length(GrangeObjectsList)){
  seqlevelsStyle(GrangeObjectsList[[i]])<- "UCSC"
}

GrangeObjectsList[[1]][1:2,]

saveRDS(GrangeObjectsList, "/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/GrangeObjectsList.rds")


# The annotatePeak function accepts a GRanges object of the regions to annotate,
# a TXDB object for gene locations and a database object name to retrieve gene names from.

AnnotatedPeakList<-list()

for(i in 1:length(GrangeObjectsList)){
  AnnotatedPeakList[[i]] <- annotatePeak(GrangeObjectsList[[i]], tssRegion = c(-2000, 200), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
    annoDb = "org.Mm.eg.db")
  }

AnnotatedPeakList

saveRDS(AnnotatedPeakList,"/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/AnnotatedPeakList.rds")
AnnotatedPeakList<-readRDS("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/AnnotatedPeakList.rds")


#=================================================================================

#Visualising peak annotation

#=================================================================================

plotAnnoBar(peakAnno)
# Plot the distribution around the tss
plotDistToTSS(peakAnno)
# Plot describing overlap between annotations

fileNames <-list.files()
SampleNames <-gsub("_peaks.xls","",fileNames)


#For individual graphs
for(i in 1:nrow(SRRList)){
  pdf(paste0(SRRList[i,],"Peak_annotation.pdf"), w = 10, h=5)
  print(plotAnnoBar(AnnotatedPeakList[i]))
  print(plotDistToTSS(AnnotatedPeakList[i]))
  dev.off()
}

names(AnnotatedPeakList)<-SampleNames

pdf(("/rds/general/user/ah3918/ephemeral/MARREFRESH/GRANGES/Graphs/Peak_annotation.pdf"), w = 10, h=5)
plotAnnoBar(AnnotatedPeakList)
plotDistToTSS(AnnotatedPeakList)
dev.off()


PeakAnnoGENESList=list()
for(i in 1:length(AnnotatedPeakList)){
  PeakAnnoGENESList[[i]]=as.data.frame(AnnotatedPeakList[[i]]);head(PeakAnnoGENESList[i])
}


length(unique(PeakAnnoGENESList[[1]]$SYMBOL))

PeakAnnoGENESList[[1]][PeakAnnoGENESList[[1]]$distanceToTSS < 200,]

for(i in 1:length(PeakAnnoGENESList)){
  PeakAnnoGENESList[[i]]=PeakAnnoGENESList[[i]][PeakAnnoGENESList[[i]]$distanceToTSS < 200,]
  PeakAnnoGENESList[[i]]=PeakAnnoGENESList[[i]][PeakAnnoGENESList[[i]]$distanceToTSS> -2000,]
}

for(i in 1:length(PeakAnnoGENESList)){
  pdf("/rds/general/user/ah3918/ephemeral/MARREFRESH/GRANGES/Graphs/Plot_distance_to_TSS.pdf")
  hist(PeakAnnoGENESList[[i]]$distanceToTSS, breaks = 100)
	dev.off()
}

head(PeakAnnoGENESList[[1]][PeakAnnoGENESList[[1]]$distanceToTSS == 0,]) # 3191

length(unique(PeakAnnoGENESList[[1]]$SYMBOL)) # 4008



#Sort based on fold_enrichment
for(i in 1:length(PeakAnnoGENESList)){
  PeakAnnoGENESList[[i]]<-PeakAnnoGENESList[[i]][order(PeakAnnoGENESList[[i]]$fold_enrichment, decreasing=TRUE),]
}



# Target genes

TargetGenes<-list()

for(i in 1:length(PeakAnnoGENESList)){
  TargetGenes[[i]]<-unique(PeakAnnoGENESList[[i]]$SYMBOL)
}
for(i in 1:length(TargetGenes)){
  TargetGenes[[i]]<-as.data.frame(TargetGenes[[i]])
}



Top500<-list()
for(i in 1:length(TargetGenes)){
  Top500[[i]]<-TargetGenes[[i]][1:500]
}


names(Top500)<-SampleNames

fileNames <-list.files()
SampleNames <-gsub("_peaks.xls","",fileNames)
names(TargetGenes)<-SampleNames





for(i in 1:length(TargetGenes)){
  write.table(TargetGenes[[i]], file = paste0("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/TargetGenes/",(SampleNames[i]),"_targets.txt",sep=""), col.names = F, row.names = F)
}
for(i in 1:length(Top500)){
  write.table(Top500[[i]], file = paste0("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/TargetGenesTop500/",(SampleNames[i]),"_top_500_targets.txt",sep=""), col.names = F, row.names = F)
}

saveRDS(TargetGenes,"/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/TargetGenes/TargetGenesRANKEDByFoldEnrichment.rds")

#=============== Collate all SRR duplicates ==================

TargetGenes<-readRDS("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/TargetGenes/TargetGenesRANKEDByFoldEnrichment.rds")

GoldStandard<-list("REST_GSE27341"=REST_GSE27341,"SOX3_Intersect"=SOX3_Intersect,
"Zic1_GSE60731"=Zic1_GSE60731,"Neurog2_GSE63620"=Neurog2_GSE63620,"Pax6Newborn_GSEGSE66961"=Pax6Newborn_GSEGSE66961,
"NeuroD2_GSE67539"=NeuroD2_GSE67539,"OTX2_GSE69724"=OTX2_GSE69724,"Olig2_GSE74646"=Olig2_GSE74646,"NKX2_1_GSE85704"=NKX2_1_GSE85704,
"LHX6_GSE85704"=LHX6_GSE85704,"SOX2_GSE90561"=SOX2_GSE90561,"SOX9_GSE117997"=SOX9_GSE117997,"Tbr2_GSE119362"=Tbr2_GSE119362)

RESTnpc_GSE27341<-Reduce(intersect,list(TargetGenes[[1]],TargetGenes[[2]],TargetGenes[[3]]))
RESTembryo_GSE27341<-Reduce(intersect,list(TargetGenes[[4]],TargetGenes[[5]],TargetGenes[[6]]))
SOX3_GSE57186<-Reduce(intersect,list(TargetGenes[[7]],TargetGenes[[8]]))
Zic1_GSE60731<-Reduce(intersect,list(TargetGenes[[9]],TargetGenes[[10]],TargetGenes[[11]]))
Neurog2_GSE63620<-TargetGenes[[15]]
Pax6Newborn_GSEGSE66961<-TargetGenes[[16]]
NeuroD2_GSE67539<-Reduce(intersect,list(TargetGenes[[17]],TargetGenes[[18]],TargetGenes[[19]]))
OTX2_GSE69724<-Reduce(intersect,list(TargetGenes[[20]],TargetGenes[[21]],TargetGenes[[22]]))
Olig2_GSE74646<-TargetGenes[[24]]
SOX3_GSE33059<-Reduce(intersect,list(TargetGenes[[27]],TargetGenes[[28]],TargetGenes[[29]]))
NKX2_1_GSE85704<-Reduce(intersect,list(TargetGenes[[30]],TargetGenes[[31]]))
LHX6_GSE85704<-TargetGenes[[32]]
SOX2_GSE90561<-TargetGenes[[33]]
SOX9_GSE117997<-Reduce(intersect,list(TargetGenes[[38]],TargetGenes[[39]],TargetGenes[[40]]))
Tbr2_GSE119362<-Reduce(intersect,list(TargetGenes[[41]],TargetGenes[[42]]))

#Remove na.action attributes from na.omit
attr(GoldStandard[1]$REST_GSE27341,"na.action")<-NULL
attr(GoldStandard[2]$SOX3_Intersect,"na.action")<-NULL

saveRDS(GoldStandard,"/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/GoldStandard.rds")

#=========== REMOVED BASED ON HOMER =================================
Tbr2_GSE63620<-TargetGenes[[14]]
Myt1lb_GSE72121<-TargetGenes[[23]]
GMNN_GSE77246<-TargetGenes[[25]]
Zic1_GSE77246<-TargetGenes[[26]]
HIF1a_GSE106220<-Reduce(intersect,list(TargetGenes[[34]])
SOX3_GSE117997<-Reduce(intersect,list(TargetGenes[[36]],TargetGenes[[37]]))




for(i in 1:length(TFList)){
  write.table(TFList[[i]], file = paste0("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/TargetGenes/",(names(TFList[i])),"_targets.txt",sep=""), col.names = F, row.names = F)
}

write.table(SampleNames,file="SRRList.txt",col.names = F, row.names = F)


TargetGenes<-readRDS("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/MARREFRESH/GRANGES/TargetGenes/TargetGenesRANKEDByFoldEnrichment.rds")


NumberUniqueGenes<-list()

test<-data.frame()

for(i in 1:length(TargetGenes)){
  print(names(TargetGenes[i]))
  print(length(TargetGenes[[i]]))
  test[i,1]<-c(length(TargetGenes[[i]]))
  test[i,2]<-c(names(TargetGenes[i]))
}

write.table(test, "Number_of_unique_reads.txt")


#Intersectsox,

SOX3_GSE57186<-Reduce(intersect,list(TargetGenes[[7]],TargetGenes[[8]]))
SOX3_GSE33059<-Reduce(intersect,list(TargetGenes[[27]],TargetGenes[[28]],TargetGenes[[29]]))

TargetList = list(a = TargetGenes[[7]], b = TargetGenes[[8]], c = TargetGenes[[27]], d = TargetGenes[[28]],
  e = TargetGenes[[29]])

str(TargetList)
TargetGenes[[7]]

AllTargets = union(TargetList[[7]],TargetList[[8]])

for(i in 3:length(TargetList)){
  AllTargets=union(AllTargets,TargetList[[i]])
}
length(AllTargets)

TargetMatrix  = matrix(nrow = length(AllTargets), ncol = 5); rownames(TargetMatrix)=AllTargets; colnames(TargetMatrix) = names(TargetList)
head(TargetMatrix)

for (i in 1:length(AllTargets)){
  for (j in 1:length(TargetList)){

    TargetMatrix[i,j] = any(AllTargets[i]==TargetList[[j]])

  }
}
head(TargetMatrix)

SummaryTarget = apply(TargetMatrix,1, sum)

SummaryTarget = SummaryTarget[SummaryTarget >= 3]; head(SummaryTarget); length(SummaryTarget)
SummaryTarget = names(SummaryTarget); head(SummaryTarget)

SummaryTarget<-na.omit(SummaryTarget)

SOX3_Intersect<-SummaryTarget









#Fold enrichment is 26.24

#checking that after ordering its still the same

#1. Liisi will send me the table enviseaged
#2. (find a reference to what fold enrichment)
#3. Create a List for each object with eeach gene ranked per fold enrichment
#4. After generating table
#5. In final product, take average fold enrichment, and keep all genes

#6. run Julia on other cell types
#7. run Liisis RScript
