options(stringsAsFactors=F)
library(colorout)
source("~/Dropbox/tools/R_functions/Liisi_functions.R")
source("~/Dropbox/tools/R_functions/000.R_functions.R")

setwd("~/GOLD_STANDARD")

Load("~PreProcessed/GS_permut_enrich.RData")
#[1,] "ResGS"
str(ResGS)


Load("~/TRRUST_permut_enrich.RData")
#[1,] "ResTRRUST"


plotEnrichPermutCell <- function(TF = "Neurog2_GSE63620", ResList = ResGS){

  par(mfrow = c(2,3))

  p = ResList$pyramidal
  p.obs = p$obs[[TF]]; names(p.obs) = "obs"
  p.diff = p$diff[[TF]]; names(p.diff) = c(1:1000); p.diff = c(p.obs, p.diff);
  p.diff = scale(p.diff)
  p.obs = p.diff[1,]; p.diff = p.diff[-1,];
  hist(p.diff, col = "#A88F07", xlab = '', xlim = c(min(p.diff),p.obs+0.5), main = "Pyramidal neurons", breaks = 15); abline(v= p.obs, lty = 5)

  i = ResList$interneurons
  i.obs = i$obs[[TF]]; names(i.obs) = "obs"
  i.diff = i$diff[[TF]]; names(i.diff) = c(1:1000); i.diff = c(i.obs, i.diff);
  i.diff = scale(i.diff)
  i.obs = i.diff[1,]; i.diff = i.diff[-1,];
  hist(i.diff, col = "#5086FF", xlab = '', xlim = c(min(i.diff),max(i.diff)+0.5), main = "Interneurons", breaks = 15); abline(v= i.obs, lty = 5)

  a = ResList$astrocytes
  a.obs = a$obs[[TF]]; names(a.obs) = "obs"
  a.diff = a$diff[[TF]]; names(a.diff) = c(1:1000); a.diff = c(a.obs, a.diff);
  a.diff = scale(a.diff)
  a.obs = a.diff[1,]; a.diff = a.diff[-1,];
  hist(a.diff, col = "#F35E5A", xlab = '', xlim = c(min(a.diff),a.obs+0.5), main = "Asotrocytes", breaks = 20); abline(v= a.obs, lty = 5)

  o = ResList$oligodendrocytes
  o.obs = o$obs[[TF]]; names(o.obs) = "obs"
  o.diff = o$diff[[TF]]; names(o.diff) = c(1:1000); o.diff = c(o.obs, o.diff);
  o.diff = scale(o.diff)
  o.obs = o.diff[1,]; o.diff = o.diff[-1,];
  hist(o.diff, col = "#1BB12B", xlab = '', xlim = c(min(o.diff),max(o.diff)+0.5), main = "Oligodendrocytes", breaks = 15); abline(v= o.obs, lty = 5)

  m = ResList$microglia
  m.obs = m$obs[[TF]]; names(m.obs) = "obs"
  m.diff = m$diff[[TF]]; names(m.diff) = c(1:1000); m.diff = c(m.obs, m.diff);
  m.diff = scale(m.diff)
  m.obs = m.diff[1,]; m.diff = m.diff[-1,];
  hist(m.diff, col = "#1EB3B7", xlab = '', xlim = c(min(m.diff),max(m.diff)), main = "Microglia", breaks = 15); abline(v= m.obs, lty = 5)

}

pdf("NeuroD2_enrichment.pdf", h = 5)
plotEnrichPermutCell(TF = "NeuroD2_GSE67539", ResList = ResGS)
dev.off()


pdf("GS_v_TRRUST_enrichment.pdf", h = 5)
par(mfrow = c(2,3))
#plotEnrichPermutTF <-function(CellType = "", ResList = ResGS, TF=c("SOX9_GSE117997","REST_GSE27341","Pax6Newborn_GSEGSE66961")){

ResList = ResGS
TF=c("SOX9_GSE117997","REST_GSE27341","Pax6Newborn_GSEGSE66961")
TFnames = c("Sox9","Rest","Pax6")

for(i in 1:length(TF)){
  p = ResList$pyramidal
  p.obs = p$obs[[TF[i]]]; names(p.obs) = "obs"
  p.diff = p$diff[[TF[i]]]; names(p.diff) = c(1:1000); p.diff = c(p.obs, p.diff);
  p.diff = scale(p.diff)
  p.obs = p.diff[1,]; p.diff = p.diff[-1,];
  hist(p.diff, col = "#A88F07", xlab = '', xlim = c(min(p.diff),p.obs+0.5), main = TFnames[i], breaks = 15); abline(v= p.obs, lty = 5)
}

ResList = ResTRRUST
TF=c("Sox9","Rest","Pax6")

for(i in 1:length(TF)){
  p = ResList$pyramidal
  p.obs = p$obs[[TF[i]]]; names(p.obs) = "obs"
  p.diff = p$diff[[TF[i]]]; names(p.diff) = c(1:1000); p.diff = c(p.obs, p.diff);
  p.diff = scale(p.diff)
  p.obs = p.diff[1,]; p.diff = p.diff[-1,];
  hist(p.diff, col = "#A88F07", xlab = '', xlim = c(min(p.diff),max(p.diff)+0.5), main = TF[i], breaks = 15); abline(v= p.obs, lty = 5)
}
dev.off()
