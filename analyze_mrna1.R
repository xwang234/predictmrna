#! /usr/bin/env Rscript

#order data by gene location, rownames of data: genesymbol
orderdatabygenelocation=function(data)
{
  genes=rownames(data)
  #Knowngenes with gene symbols downloaded from https://genome.ucsc.edu/cgi-bin/hgTables
  allgenes=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/knownGene1.txt",sep="\t")
  allgenes=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/knownGenehg18.txt",sep="\t")
  allgenes=allgenes[,c(1,2,4,14)]
  colnames(allgenes)=c("kgid","chr","start","symbol")
  allgenes$chr=as.character(allgenes$chr)
  chr24=paste0("chr",1:24)
  idx=allgenes$chr=="chrX"
  allgenes$chr[idx]="chr23"
  idx=allgenes$chr=="chrY"
  allgenes$chr[idx]="chr24"
  idx=allgenes$chr %in% chr24
  allgenes=allgenes[idx,]
  allgenes$chr=gsub("chr","",allgenes$chr) 
  #order allgenes
  allgenes1=NULL
  for (i in 1:24)
  {
    tmp=allgenes[allgenes$chr==i,]
    idx=order(tmp$start)
    allgenes1=rbind(allgenes1,tmp[idx,])
  }
  
  allgenes1=cbind(allgenes1,1:nrow(allgenes1))
  colnames(allgenes1)[ncol(allgenes1)]="allgenesID"
  idx=genes %in% allgenes1[,4] #10630 genes were ordered
  #genes match symbols in Knowngenes table
  genes1=genes[idx]
  data1=data[genes1,]
  #genes not match, they may have alias gene name
  genes1table=data.frame(matrix(nrow=length(genes1),ncol=2))
  colnames(genes1table)=c("symbol","ID")
  genes1table[,1]=genes1
  genes1table[,2]=1:length(genes1)
  tmp=merge(allgenes1,genes1table,by="symbol")
  tmp=tmp[order(tmp[,"allgenesID"]),] #ordered by Knowngenes table
  ordID=unique(tmp[,"ID"]) #Use ID in genes1
  genes2=genes[!idx] #1233 genes were not ordered, add to the list
  data2=data[genes2,]
  res=rbind(data1[ordID,],data2)
}

#get the chromosome of genes
getgenechr=function(genes)
{
  chr24=paste0("chr",1:24)
  chrs=rep(25,length(genes)) #unknown ones set as 25
  allgenes=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/knownGene1.txt",sep="\t")
  allgenes=allgenes[,c(1,2,4,14)]
  colnames(allgenes)=c("kgid","chr","start","symbol")
  allgenes$chr=as.character(allgenes$chr)
  idx=allgenes$chr=="chrX"
  allgenes$chr[idx]="chr23"
  idx=allgenes$chr=="chrY"
  allgenes$chr[idx]="chr24"
  chrs=sapply(genes, function(gene) {
    idx=allgenes$symbol==gene
    idx1=allgenes$chr[idx] %in% chr24
    if (sum(idx1)>0)
    {
      res=allgenes$chr[idx][idx1]
      res=res[1]
      res=as.numeric(substr(res,4,nchar(res)))
    }else
    {
      res=25
    }
    return(res)
  })
  return(chrs)
}

gengenecolor=function(genes)
{
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  genechr=getgenechr(genes)
  res=color2[genechr]
}
#replace NA to 0
NA2zero=function(data)
{
  res=data
  if (dim(data)[1]<dim(data)[2])
  {
    for (i in 1:nrow(data))
    {
      idx=is.na(data[i,])
      res[i,idx]=0
    }
  }else
  {
    for (i in 1:ncol(data))
    {
      idx=is.na(data[,i])
      res[idx,i]=0
    }
  }
  
  return(res)
}

#draw the correlation of mrna and methylation
load("/fh/fast/dai_j/CancerGenomics/Ovarian/mrna_copynumber_methylation_mutation.RData")
mrna=orderdatabygenelocation(mrna)
methylation=orderdatabygenelocation(methylation)
methylation=fillNAs(methylation)
methylation1=orderdatabygenelocation(methylation1)
methylation1=fillNAs(methylation1)
copynumber=orderdatabygenelocation(copynumber)
copynumber=fillNAs(copynumber)
mutation=orderdatabygenelocation(mutation)
#compute correlation of x with all rows of Y
corxy=function(x,y)
{
  x=as.numeric(x)
  y=as.numeric(y)
  res=format(cor(x,y),digits=2)
}
plinecor=function(x,Y)
{
  x=as.numeric(x)
  Y=as.matrix(Y)
  res=rep(NA,nrow(Y))
  res=apply(Y,1,corxy,y=x)
}

#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
njob=100
library(Rmpi)
#mpi.spawn.Rslaves(nslaves=4)
#mpi.spawn.Rslaves()
mpi.spawn.Rslaves(needlog = FALSE)
# In case R exits unexpectedly, have it automatically clean up # resources taken up by Rmpi (slaves, memory, etc...)
.Last <- function()
{ if (is.loaded("mpi_initialize")){
  if (mpi.comm.size(1) > 0){
    print("Please use mpi.close.Rslaves() to close slaves.")
    mpi.close.Rslaves()
  }
  print("Please use mpi.quit() to quit R")
  .Call("mpi_finalize")
}
}
mpi.bcast.Robj2slave(mrna)
mpi.bcast.Robj2slave(methylation)
mpi.bcast.Robj2slave(methylation1)
mpi.bcast.Robj2slave(mutation)
mpi.bcast.Robj2slave(copynumber)
mpi.bcast.Robj2slave(corxy)
mpi.bcast.Robj2slave(plinecor)

#cor(data1,data2)
cormatrix=function(data1,data2,output)
{
  res1 <- NULL
  nrun <- ceiling(nrow(data1)/1000)
  for (j in 1:nrun){
    cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):nrow(data1)
    res=mpi.parRapply(X=data1[cseq,],FUN=plinecor,Y=data2,job.num=njob)
    tmp=matrix(unlist(res),byrow=TRUE,ncol=nrow(data1))
    res1=rbind(res1,tmp)
    rownames(res1)[cseq]=rownames(data1)[cseq]
    colnames(res1)=rownames(data2)
    write.table(res1,file=output,col.names=T,row.names=T,sep="\t",quote=F)
  }
  return(res1)
}

data1=mrna
data2=copynumber
output="corr_mrna_copynumber.txt"
res=cormatrix(data1,data2,output)

data1=mrna
data2=methylation
output="corr_mrna_methylation.txt"
res1=cormatrix(data1,data2,output)

data1=mrna
data2=methylation1
output="corr_mrna_methylation1.txt"
res2=cormatrix(data1,data2,output)

data1=mrna
data2=mutation
output="corr_mrna_mutation.txt"
res3=cormatrix(data1,data2,output)

#mpi.close.Rslaves()
mpi.close.Rslaves(dellog = FALSE)
mpi.quit()

#draw heatmap
cor_mrna_copynumber=read.table(file="corr_mrna_copynumber.txt",header=T,sep="\t")
cor_mrna_methylation=read.table(file="corr_mrna_methylation.txt",header=T,sep="\t")
cor_mrna_methylationn=read.table(file="corr_mrna_methylation1.txt",header=T,sep="\t")
cor_mrna_mutation=read.table(file="corr_mrna_mutation.txt",header=T,sep="\t")

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

genes=rownames(mrna)
genecolor=gengenecolor(genes)
printheatmap=function(cordata,output,qt)
{
  #my_palette <- colorRampPalette(c("green","yellow","red"))(n = 299)
  #used for colorpalette, 3cor1,6cor2,3cor3. work with breaks, so that breaks divided into 3 segments, each segment has 3 subsegments.
  #in the first segment, color scheme is col1-col1,col1-col1,col1-col2, in the thrid subsegment color finaly sets to col2.
  #in the second segment, color is set as col2, col2-col2,col2-col2,col2-col2.
  formcolors=function(colors)
  {
    allcolors=c(rep(colors[1],3),rep(colors[2],6),rep(colors[3],3))
  }
  my_palette <- colorRampPalette(formcolors(c("green","white","red")))(n = 299)
  #qt=c(0,0.01,0.99,1)
  #qt=c(0,0.0001,0.01,0.99,0.9999,1) #show extreme values
  tmp=as.matrix(quantile(cordata,qt,na.rm=T))
#   sec1=seq(tmp[[1]]-0.01,tmp[[2]],length=100)
#   sec2=seq(tmp[[2]]+mean(diff(sec1)),tmp[[length(qt)-1]]-mean(diff(sec1)),length=100)
#   sec3=seq(tmp[[length(qt)-1]],tmp[[length(qt)]]+0.01,length=100)
  if (tmp[[2]]==tmp[[length(qt)-1]])
  {
    tmp[[2]]=tmp[[2]]-0.01
  }
  sec1=seq(tmp[[1]],tmp[[2]],length=100)
  sec2=seq(tmp[[2]]+0.0001,tmp[[length(qt)-1]],length=100)
  sec3=seq(tmp[[length(qt)-1]]+0.0001,tmp[[length(qt)]],length=100)
  col_breaks=c(sec1,sec2,sec3)
  #png(filename = output,width = 500, height = 500, units = "px")
  png(filename = output,width = 800, height = 800, units = "px")
  heatmap.2(as.matrix(cordata),
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            #margins =c(4,2),     # widens margins around plot
            breaks=col_breaks,    # enable color transition at specified limits
            symkey=F,             #Boolean indicating whether the color key should be made symmetric about 0
            #keysize=2,            #numeric value indicating the size of the key
            col=my_palette,       # use on color palette defined earlier 
            dendrogram="none",
            Rowv=NULL,
            Colv=NULL,
            labCol="",
            cexCol=1,
            labRow="",
            ColSideColors=genecolor,
            RowSideColors=genecolor,
            cex.axis=1,
            xlab="",ylab=""
  )
  dev.off()
}
qt=c(0,0.01,0.99,1)
qt=c(0,0.1,0.9,1)
qth=c(0,0.001,0.01,0.99,0.999,1) #show extreme values, could use for lasso cases

cor_mrna_copynumber1=NA2zero(cor_mrna_copynumber)
cordata=cor_mrna_copynumber1
output="cor_mrna_copynumber.png"
printheatmap(cordata,output,qt)
#output="cor_mrna_copynumberh.png"
#printheatmap(cordata,output,qth)

cor_mrna_methylation1=NA2zero(cor_mrna_methylation)
cordata=cor_mrna_methylation1
output="cor_mrna_methylation.png"
printheatmap(cordata,output,qt)
#output="cor_mrna_methylationh.png"
#printheatmap(cordata,output,qth)

cor_mrna_methylationn1=NA2qtho(cor_mrna_methylationn)
cordata=cor_mrna_methylationn1
output="cor_mrna_methylationn.png"
printheatmap(cordata,output,qt)
#output="cor_mrna_methylationnh.png"
#printheatmap(cordata,output,qth)

cor_mrna_mutation1=NA2zero(cor_mrna_mutation)
cordata=cor_mrna_mutation1
output="cor_mrna_mutation.png"
printheatmap(cordata,output,qt)
#output="cor_mrna_mutationh.png"
#printheatmap(cordata,output,qth)


#cross anaylize predictmrna and correlation result
which(rownames(predictmrna1se)=="NOC2L")
#[1] 6976
tmp=as.character(predictmrna1se[6976,492])
#1] "781,5913,12684,13431,13920,15466,16711,20085,21489,22712"
tmp=as.integer(unlist(strsplit(tmp,',')))
colnames(X)[tmp]
#[1] "CP_B3GALT2" "CP_NOC2L"   "Me_EDN2"    "Me_GH2"     "Me_HMGCL"   "Me_MYH1"    "Me_PRKRIR"  "MU_ARL6IP5" "MU_FILIP1L" "MU_MED12"
count1se[6976,]
#      CP ME MU
#NOC2L  2  5  3
plot(as.numeric(cor_mrna_copynumber1[1,]))

#calculate correlation using LASSO

#fill missing value with row mean
fillNAs=function(data)
{
  #rowmean=rowMeans(data,na.rm=TRUE)
  rowmean=rowSums(data,na.rm=T)/ncol(data)
  for (i in 1:nrow(data))
  {
    idx=is.na(data[i,])
    if (sum(idx)>0)
    {
      data[i,idx]=rowmean[i]
    }
  }
  return(data)
}

cv.glmnets=function(X,Y,nfolds=10,ntime=10)
{
  res=data.frame(matrix(NA,ncol=2,nrow=1))
  colnames(res)=c('lambda.min','lambda.1se')
  
  #ntime:number of times run cv.glmnet, glmnet is not stable when p is large and nfolds is small. using leave one out is time-consuming
  lambda=NULL
  lambdamin=NULL
  lambda1se=NULL
  for (i in 1:ntime){
    cv <- cv.glmnet(X, Y,nfolds=nfolds)
    if (length(lambda)<length(cv$lambda)) lambda=cv$lambda #record the longest lambda, length of lambda varies in 96-99
    lambdamin=c(lambdamin,cv$lambda.min)
    lambda1se=c(lambda1se,cv$lambda.1se)
  }
  temp=as.matrix(table(lambdamin))
  temp=temp[order(temp[,1],decreasing = TRUE),1]
  maxtime=temp[1]
  idx=which(temp==maxtime) #the most frequently appeared lambdas
  bestlambdas=as.numeric(names(temp[idx]))
  res$lambda.min=max(bestlambdas)
  
  temp=as.matrix(table(lambda1se))
  temp=temp[order(temp[,1],decreasing = TRUE),1]
  maxtime=temp[1]
  idx=which(temp==maxtime) #the most frequently appeared lambdas
  bestlambdas=as.numeric(names(temp[idx]))
  res$lambda.1se=max(bestlambdas)
  
  return(res)
}

lassoplinecor=function(linedata1,data2) #data2, data1:ncol,number of samples
{
  res1=data.frame(matrix(NA,nrow=1,ncol=nrow(data2)))
  rownames(res1)=rownames(linedata1)
  colnames(res1)=rownames(data2)
  res2=res1
  linedata1=as.numeric(linedata1)
  data2=t(as.matrix(data2))
  fit = glmnet(data2, linedata1)
  #cvfit=cv.glmnet(data2,linedata1)
  lambdas =cv.glmnets(data2,linedata1,nfolds=10,ntime=10)
  sel_set <- which(as.matrix(coef(fit,s=lambdas$lambda.min))[,1]!=0)
  if (length(sel_set)>1) # at least one variable selected
  {
    sel_variables=(sel_set-1)[2:length(sel_set)]
    data2sel= data2[,(sel_set-1)[2:length(sel_set)]]
    res1[sel_variables]=plinecor(linedata1,t(data2sel))
  }
  
  sel_set <- which(as.matrix(coef(fit,s=lambdas$lambda.1se))[,1]!=0)

    if (length(sel_set)>1) # at least one variable selected
  {
    sel_variables=(sel_set-1)[2:length(sel_set)]
    data2sel= data2[,(sel_set-1)[2:length(sel_set)]]
    res2[sel_variables]=plinecor(linedata1,t(data2sel))
  }
  result=list(resmin=res1,res1se=res2)
}

load("/fh/fast/dai_j/CancerGenomics/Ovarian/mrna_copynumber_methylation_mutation.RData")
mrna=orderdatabygenelocation(mrna)
methylation=orderdatabygenelocation(methylation)
methylation=fillNAs(methylation)
methylation1=orderdatabygenelocation(methylation1)
methylation1=fillNAs(methylation1)
copynumber=orderdatabygenelocation(copynumber)
copynumber=fillNAs(copynumber)
mutation=orderdatabygenelocation(mutation)
#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
njob=100
library(Rmpi)
#mpi.spawn.Rslaves(nslaves=4)
#mpi.spawn.Rslaves()
mpi.spawn.Rslaves(needlog = FALSE)
# In case R exits unexpectedly, have it automatically clean up # resources taken up by Rmpi (slaves, memory, etc...)
.Last <- function()
{ if (is.loaded("mpi_initialize")){
  if (mpi.comm.size(1) > 0){
    print("Please use mpi.close.Rslaves() to close slaves.")
    mpi.close.Rslaves()
  }
  print("Please use mpi.quit() to quit R")
  .Call("mpi_finalize")
}
}
mpi.bcast.Robj2slave(mrna)
mpi.bcast.Robj2slave(methylation)
mpi.bcast.Robj2slave(methylation1)
mpi.bcast.Robj2slave(mutation)
mpi.bcast.Robj2slave(copynumber)
mpi.bcast.Robj2slave(corxy)
mpi.bcast.Robj2slave(plinecor)
mpi.bcast.Robj2slave(lassoplinecor)
mpi.bcast.cmd(library(glmnet))
mpi.bcast.cmd(require(methods))
mpi.bcast.Robj2slave(cv.glmnets)

#cor(data1,data2)
lassocormatrix=function(data1,data2,output1,output2)
{
  res1=res2=NULL
  nrun <- ceiling(nrow(data1)/1000)
  for (j in 1:nrun){
    cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):nrow(data1)
    res=mpi.parRapply(X=data1[cseq,],FUN=lassoplinecor,data2=data2,job.num=njob)
    tmp=matrix(unlist(res),byrow=TRUE,ncol=nrow(data2)) #for each line in mrna, return 2 lines resmin and res1se
    res1=rbind(res1,tmp[seq(1,nrow(tmp),2),]) #odd row resmin
    res2=rbind(res2,tmp[seq(2,nrow(tmp),2),]) #even row res1se
    rownames(res1)[cseq]=rownames(res2)[cseq]=names(res)
    colnames(res1)=colnames(res2)=rownames(data2)
    write.table(res1,file=output1,col.names=T,row.names=T,sep="\t",quote=F)
    write.table(res2,file=output2,col.names=T,row.names=T,sep="\t",quote=F)
  }
  return(rbind(res1,res2))
}

data1=mrna
data2=copynumber
output1="corr_mrna_copynumber_lassomin.txt"
output2="corr_mrna_copynumber_lasso1se.txt"
res=lassocormatrix(data1,data2,output1,output2)

data1=mrna
data2=methylation
output1="corr_mrna_methylation_lassomin.txt"
output2="corr_mrna_methylation_lasso1se.txt"
res1=lassocormatrix(data1,data2,output1,output2)

data1=mrna
data2=methylation1
output1="corr_mrna_methylation1_lassomin.txt"
output2="corr_mrna_methylation1_lasso1se.txt"
res2=lassocormatrix(data1,data2,output1,output2)

data1=mrna
data2=mutation
output1="corr_mrna_mutation_lassomin.txt"
output2="corr_mrna_mutation_lasso1se.txt"
res3=lassocormatrix(data1,data2,output1,output2)

#mpi.close.Rslaves()
mpi.close.Rslaves(dellog = FALSE)
mpi.quit()

cor_mrna_copynumber_lasso=read.table(file="corr_mrna_copynumber_lasso1se.txt",header=T,sep="\t")
cor_mrna_methylation_lasso=read.table(file="corr_mrna_methylation_lasso1se.txt",header=T,sep="\t")
cor_mrna_methylationn_lasso=read.table(file="corr_mrna_methylation1.txt",header=T,sep="\t")
cor_mrna_mutation_lasso=read.table(file="corr_mrna_mutation_lasso1se.txt",header=T,sep="\t")

cor_mrna_copynumber_lasso1=NA2zero(cor_mrna_copynumber_lasso)
cordata=cor_mrna_copynumber_lasso1
output="cor_mrna_copynumber_lasso.png"
printheatmap(cordata,output,qt)
#output="cor_mrna_copynumberh.png"
#printheatmap(cordata,output,qth)

cor_mrna_methylation_lasso1=NA2zero(cor_mrna_methylation_lasso)
cordata=cor_mrna_methylation_lasso1
output="cor_mrna_methylation_lasso.png"
printheatmap(cordata,output,qt)
#output="cor_mrna_methylationh.png"
#printheatmap(cordata,output,qth)

cor_mrna_methylationn_lasso1=NA2qtho(cor_mrna_methylationn_lasso)
cordata=cor_mrna_methylationn_lasso1
output="cor_mrna_methylationn_lasso.png"
printheatmap(cordata,output,qt)
#output="cor_mrna_methylationnh.png"
#printheatmap(cordata,output,qth)

cor_mrna_mutation_lasso1=NA2zero(cor_mrna_mutation_lasso)
cordata=cor_mrna_mutation_lasso1
output="cor_mrna_mutation_lasso.png"
printheatmap(cordata,output,qt)
#output="cor_mrna_mutationh.png"
#printheatmap(cordata,output,qth)
