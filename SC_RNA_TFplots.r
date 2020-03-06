

# !!!!!!!!!!!!!
# Islet1
# Tbr2 cant be found!!!!
# Neither can Sox3!!!
#Create violin plots for TFs in each cell types

TFlist=read.delim("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/TFlist.txt",header=FALSE)

pdf("Gmnn_ViolinPlot.pdf")
VlnPlot(HS.data.zeisel,"Gmnn")
dev.off()



readRDS("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/SCANALYSIS/HS100bp/HS_FinalRawFiltered_dgCMatrix.rds")

x=HS.FinalRawFilteredMatrix["Pax6",]
table(x>0)


for (i in 1:nrow(TFlist)){
  pdf(paste0(TFlist[i,],"_ViolinPlot.pdf"))
  print(VlnPlot(HS.data.zeisel,TFlist[i,]))
  dev.off()
}



for (i in 1:nrow(TFlist)){
  b=HS.FinalRawFilteredMatrix[i,]
  print(table(b>0))
  print(as.character(TFlist[i,]))
}

#Create a data frame for it
FinalDataFrame<-as.data.frame(matrix(0,ncol=20,nrow=2))


for (i in 1:nrow(TFlist)){
  b=HS.FinalRawFilteredMatrix[i,]
  FinalDataFrame[i]<-cbind(table(b>0))
}


FinalDataFrame<-t(FinalDataFrame)

for (c in 1:nrow(TFlist)){
  row.names(FinalDataFrame)[c]<-as.matrix(TFlist[c,])
}

colnames(FinalDataFrame)<-c("True.Zeros","Expression.Vals>0")

XataFrame<-cbind(rownames(FinalDataFrame),FinalDataFrame[,1],FinalDataFrame[,2])
FinalDataFrame2<-as.data.frame(XataFrame)
colnames(FinalDataFrame2)<-c("TF","True.Zeros","Expression.Vals.Above.Zero")

pdf("Testggplot.pdf")
ggplot(testmelt, aes(TF, Freq))+
   geom_histogram(position="dodge")
   dev.off()

dfz<-melt(FinalDataFrame2[,c("True.Zeros","Expression.Vals.Above.Zero")],id.vars=1)

ggplot(dfz,aes(x = ,y = value)) +
    geom_bar(aes(fill = variable),stat = "identity",position = "dodge") +
    scale_y_log10()

          Input Rtime Rcost Rsolutions  Btime Bcost
   1   12-proc.     1    36     614425     40    36
   2   15-proc.     1    51     534037     50    51
   3    18-proc     5    62    1843820     66    66
   4    20-proc     4    68    1645581 104400    73
   5 20-proc(l)     4    64    1658509  14400    65
   6    21-proc    10    78    3923623 453600    82,header = TRUE,sep = "")

   dfm <- melt(df[,c('Input','Rtime','Btime')],id.vars = 1)

   ggplot(dfm,aes(x = Input,y = value)) +
       geom_bar(aes(fill = variable),stat = "identity",position = "dodge") +


tbl <- with(FinalDataFrame, table(Species, Depth))
