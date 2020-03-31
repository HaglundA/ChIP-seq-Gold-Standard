options(stringsAsFactors=F)
rm(list=ls())
library(colorout)

#----------------------------------------------- Functions ----------------------------------------------------------------

# Function for viewing matrices and data frames with more information
Head<-function(data_obj,nlist=1,ncol=1:5,nrow=1:10){
  cat("\n\tobject class : ",class(data_obj),"\n\n")
  prev=""

  if(class(data_obj)=="list"){

  	if(length(names(data_obj))<50){
	    cat("\t\tlist contains",length(names(data_obj)),"objects:\n")
	    cat("\t\t\t",as.matrix(names(data_obj)),"\n\n",sep="   ")
	    cat("\t\tlist[[",nlist,"]] contains ",class((data_obj[[1]]))," :",sep="")

	    prev=paste("list[[",nlist,"]] contains :",sep="")
	    data_obj=data_obj[[nlist]]
    }

  	if(length(names(data_obj))>50){
	    cat("\t\tlist contains",length(names(data_obj)),"objects, first 50:\n")
	    cat("\t\t\t",as.matrix(names(data_obj)[1:min(50,length(names(data_obj)))]),"...\n\n",sep="   ")
	    cat("\t\tlist[[",nlist,"]] contains ",class((data_obj[[1]]))," :",sep="")

	    prev=paste("list[[",nlist,"]] contains :",sep="")
	    data_obj=data_obj[[nlist]]
    }

  }

  if(class(data_obj)=="list"){
    cat("",length(names(data_obj)),"objects:\n")
    cat("\t\t\t",as.matrix(names(data_obj)),"\n\n",sep="   ")
  }

  if(class(data_obj)=="data.frame" | class(data_obj)=="matrix"){
    cat("\n\n")
    print(data_obj[min(1,min(nrow)):min(max(nrow),nrow(data_obj)),min(1,min(ncol)):min(max(ncol),ncol(data_obj)),drop=F])
    cat("\n\t",prev,class(data_obj),"dimensions : ",dim(data_obj),"\n")
    cat("\t\tis.numeric :",is.numeric(data_obj))

    if(is.numeric(data_obj)){
    	cat('\tmin=',min(data_obj,na.rm=T),'max=',max(data_obj,na.rm=T),'\n')
    }
    cat('\n')
  }

  if(class(data_obj)=="vector"){
    cat("\n\n")
    print(data_obj[nrow])
    cat("\n\t",prev,class(data_obj),"length : ",length(data_obj),"\n")
    cat("\t\tis.numeric :",is.numeric(data_obj))

    if(is.numeric(data_obj)){
    	cat('\tmin=',min(data_obj,na.rm=T),'max=',max(data_obj,na.rm=T),'\n')
    }
    cat('\n')

  }

  if(class(data_obj)!="list" & class(data_obj)!="data.frame" & class(data_obj)!="matrix" & class(data_obj)!="vector"){
    cat("\n\n")
    str(data_obj)
    cat("\n\tis.numeric :",is.numeric(data_obj))
    if(is.numeric(data_obj)){
    	cat('\tmin=',min(data_obj,na.rm=T),'max=',max(data_obj,na.rm=T),'\n')
    }
    cat('\n')

  }
}

# custom fucntion to write tables with new defaults
write.delim<-function(mat,file,row.names=T,col.names=T,missing.value.char="NA",sep="\t",...){
  if(col.names==T & row.names==T){
    col.names=NA
  }
  write.table(mat,file,row.names=row.names,col.names=col.names,sep=sep,quote=F,na=missing.value.char,...)
}


# scale the dataset to have values between 0 and 1
range01=function(x){(x-min(x))/(max(x)-min(x))}


#>> comparing the working input data to the Astrocyte data
yeast3 = read.delim('~/50_yeast3_large.txt', row.names = 1)
yeast3 = as.matrix(yeast3)
Head(yeast3)
# 50 2100
# no column headers
# genes as rows
# values between 0 and 1

sum(yeast3[1,])
sum(yeast3[,1])
# rowsum > colsum

table(yeast3 == 0)
#FALSE   TRUE
#104744    256

# >> single cell astro data
astro = read.delim("~/Astrocyte_Raw_Expression_Matrix.txt", row.names =1)
#astro = read.delim("/Users/ll3515/Dropbox/Academic/Students/Alexander_Haglund/Project/NETWORKINFERENCE/Astrocyte_Raw_Expression_Matrix.txt", row.names = 1)
astro = as.matrix(astro)
Head(astro)
#19584 631

table(astro == 0)
#FALSE     TRUE
# 696486 11661018

# number of non-zero counts per gene
genecount = apply(astro, 1, function(x){sum(x!=0)})
hist(genecount)
table(genecount == 0)
#FALSE  TRUE
#17120  2464

#>> is the median expression lower in genes that have more drop-outs?
#astro = astro[genecount!=0,]
#genecount = apply(astro, 1, function(x){sum(x!=0)})
#Head(genecount)

# calculating mean of genes with non-zero values
#non0mean = apply(astro, 1, function(x){mean(x[x!=0])})
#Head(non0mean)

# calculating mean with all values included
#absMean = apply(astro, 1, mean)

#par(mfrow=c(3,1))
#plot(genecount,non0mean); abline(v = 63.1, col = "red")
#plot(absMean,non0mean);
#plot(genecount,absMean);


#>> creating a toy dataset with 80% of cells expressing any given gene
keep=genecount[genecount >= ncol(astro)*0.8]; length(keep)

toy = astro[names(keep),]; Head(toy)

# cells with non-zero expression
tmp = apply(toy, 2, function(x){sum(x!=0)})
any(tmp==0)
# all of the cells are expressing at least one of the genes in the gene set


for(i in 1:nrow(toy)){

  toy[i,] = range01(toy[i,])

}

Head(toy)
# 73 631
# genes as rows
# values between 0 and 1

sum(toy[1,])
sum(toy[,1])
# rowsum > colsum

setwd("/simulated_datasets")
write.delim(toy, col.names = F, file = "astroToy.txt")
# ran the toy example in julia successfully


# in a supposedly relatively homogeneous cell population, what is the proportion of cells that should express any given gene?
# based on previous gene expression based network inference methods using mutual information,
# the reccomended number of samples expressing a any given gene is 100.

keep=genecount[genecount >= 100]
#1785

astroData = astro[names(keep),]; Head(astroData)

# remove all the mitochondrial genes
mito = rownames(astroData)[grep("mt-.*",rownames(astroData))]; length(mito) # 13

astroData = astroData[!rownames(astroData)%in%mito,]; Head(astroData)

# cells with non-zero expression
tmp = apply(astroData, 2, function(x){sum(x!=0)})
any(tmp==0)
# all of the cells are expressing at least one of the genes in the gene set


# scale the dataset to have values between 0 and 1
for(i in 1:nrow(astroData)){
  astroData[i,] = range01(astroData[i,])
}
Head(astroData)
# 1772 631
# genes as rows
# values between 0 and 1

write.delim(astroData, col.names = F, file = "astroDataPreProcessed.txt")
