#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
#segid=as.integer(args[1])
njob=as.integer(args[1])
setwd("/fh/fast/dai_j/CancerGenomics/Ovarian")

load("./mrna_copynumber_methylation_mutation.RData")
library(glmnet)
require(methods)
#remove NA rows
rmNArows=function(data)
{
  
  idx=apply(data,1,function(x) {res=TRUE 
  if (sum(is.na(x))==length(x)) res=FALSE
  res})
  data=data[idx,]
}

#remove all-zero rows
rmzerorows=function(data)
{
  
  idx=apply(data,1,function(x) {res=TRUE 
  if (sum(x==0)==length(x)) res=FALSE
  res})
  data=data[idx,]
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

#fill missing value with mean
fillNAs=function(data)
{
  rowmean=rowMeans(data,na.rm=TRUE)
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
copynumber1=rmNArows(copynumber)
copynumber1=NA2zero(copynumber1)
rownames(copynumber1)=paste0("CP_",rownames(copynumber1))
methylation1=rmNArows(methylation)
methylation1=fillNAs(methylation1)
rownames(methylation1)=paste0("Me_",rownames(methylation1))
mutation1=rmzerorows(mutation)
rownames(mutation1)=paste0("MU_",rownames(mutation1))

#form X matrix
#copynumber, methylation, mutation
X=cbind(t(as.matrix(copynumber1)),t(as.matrix(methylation1)),t(as.matrix(mutation1)))
# #Check columns contain many NA
# numNA=apply(X,2,function(x) sum(is.na(x)))
# sum(numNA>dim(X)[1]/10)
# #[1] 35: 35 features have NAs in more than 10% of samples
# idx=numNA>dim(X)[1]/10
# X=X[,!idx]
# 
# #Check rows contain many NA
# numNA=apply(X,1,function(x) sum(is.na(x)))
# sum(numNA>dim(X)[2]/10)
# #[1] 1: 1 sample TCGA-61-1721 have 10388 NAs, more than 10% of features
# idx=numNA>dim(X)[2]/10
# X=X[!idx,]
# 
# numNA=apply(X,1,function(x) sum(is.na(x)))
# sum(numNA>0)
# #[1] 332 :number of samples have NA
# numNA=apply(X,2,function(x) sum(is.na(x)))
# sum(numNA>0)
# #[1] 782 :number of features have NA

library(Rmpi)
#mpi.spawn.Rslaves(nslaves=4)
mpi.spawn.Rslaves()
#mpi.spawn.Rslaves(needlog = FALSE)

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

predictrmrna=function(rowid,mrna,X)
{
  require(glmnet)
  require(methods)
  #res1=data.frame(matrix(NA,nrow=nrow(mrna),ncol=ncol(mrna)+4))
  res1=data.frame(matrix(NA,nrow=length(rowid),ncol=dim(mrna)[2]+4))
  #res1=data.frame(matrix(NA,nrow=11864,ncol=493))
  rownames(res1)=rownames(mrna)[rowid]
  colnames(res1)=c(colnames(mrna),"rsquared","numvariables","selvariables","cvwork")
  res2=res1
  
  count=1
  for (id in rowid)
  {
    Y=t(as.matrix(mrna[id,]))
    # Y=Y[!idx]
    fit = glmnet(X, Y)
    #   cvfit <- cv.glmnet(X,Y,nfolds=length(Y)) #time-consuming leave one out
    #   cvfit$lambda.1se
    #   cvfit$lambda.min
    #   sel_set <- which(as.matrix(coef(fit,s=cvfit$lambda.1se))[,1]!=0)  
    lambdas =cv.glmnets(X,Y,nfolds=10,ntime=10)
    sel_set <- which(as.matrix(coef(fit,s=lambdas$lambda.min))[,1]!=0)
    if (length(sel_set)>1) # at least one variable selected
    {
      sel_variables=paste((sel_set-1)[2:length(sel_set)],collapse="," )
      Xsel= X[,(sel_set-1)[2:length(sel_set)]]
      fit1=lm(Y~Xsel)
      #summary(fit1)
      rsquared=summary(fit1)$r.squared
      Yhat=fit1$fitted.values
      res1[count,]=c(Yhat,rsquared,length(sel_set)-1,sel_variables,1)
    }else
    {
      sel_set <- which(as.matrix(coef(fit,s=0.07))[,1]!=0) #cv.glmnet select no variables, set a specific lambda
      if (length(sel_set)>1) # at least one variable selected
      {
        sel_variables=paste((sel_set-1)[2:length(sel_set)],collapse="," )
        Xsel= X[,(sel_set-1)[2:length(sel_set)]]
        fit1=lm(Y~Xsel)
        #summary(fit1)
        rsquared=summary(fit1)$r.squared
        Yhat=fit1$fitted.values
        res1[count,]=c(Yhat,rsquared,length(sel_set)-1,sel_variables,0)
      }
    }
    
    
    sel_set <- which(as.matrix(coef(fit,s=lambdas$lambda.1se))[,1]!=0)
    if (length(sel_set)>1) # at least one variable selected
    {
      sel_variables=paste((sel_set-1)[2:length(sel_set)],collapse="," )
      Xsel= X[,(sel_set-1)[2:length(sel_set)]]
      fit1=lm(Y~Xsel)
      #summary(fit1)
      rsquared=summary(fit1)$r.squared
      Yhat=fit1$fitted.values
      res2[count,]=c(Yhat,rsquared,length(sel_set)-1,sel_variables,1)
    }else
    {
      sel_set <- which(as.matrix(coef(fit,s=0.1))[,1]!=0) #cv.glmnet select no variables, set a specific lambda
      if (length(sel_set)>1) # at least one variable selected
      {
        sel_variables=paste((sel_set-1)[2:length(sel_set)],collapse="," )
        Xsel= X[,(sel_set-1)[2:length(sel_set)]]
        fit1=lm(Y~Xsel)
        #summary(fit1)
        rsquared=summary(fit1)$r.squared
        Yhat=fit1$fitted.values
        res2[count,]=c(Yhat,rsquared,length(sel_set)-1,sel_variables,0)
      }
    }
    if (count%%50==0) print(count)
    count=count+1
  }
  result=list(resmin=res1,res1se=res2)
  return(result)
}

predictrmrna_1row=function(z,Xall)
{
  
  require(glmnet)
  require(methods)
  #res1=data.frame(matrix(NA,nrow=nrow(mrna),ncol=ncol(mrna)+4))
  
  res1=data.frame(matrix(NA,nrow=1,ncol=length(z)+4))
  #res1=data.frame(matrix(NA,nrow=11864,ncol=493))
  rownames(res1)=rownames(z)
  colnames(res1)=c(colnames(z),"rsquared","numvariables","selvariables","cvwork")
  res2=res1
  
  count=1
  #for (id in rowid)
  #{
    #Y=t(as.matrix(z))
    Y=as.numeric(z)
    # Y=Y[!idx]
    fit = glmnet(Xall, Y)
    #   cvfit <- cv.glmnet(X,Y,nfolds=length(Y)) #time-consuming leave one out
    #   cvfit$lambda.1se
    #   cvfit$lambda.min
    #   sel_set <- which(as.matrix(coef(fit,s=cvfit$lambda.1se))[,1]!=0)  
    lambdas =cv.glmnets(Xall,Y,nfolds=10,ntime=10)
    sel_set <- which(as.matrix(coef(fit,s=lambdas$lambda.min))[,1]!=0)
    if (length(sel_set)>1) # at least one variable selected
    {
      sel_variables=paste((sel_set-1)[2:length(sel_set)],collapse="," )
      Xsel= Xall[,(sel_set-1)[2:length(sel_set)]]
      fit1=lm(Y~Xsel)
      #summary(fit1)
      rsquared=summary(fit1)$r.squared
      Yhat=fit1$fitted.values
      res1[count,]=c(Yhat,rsquared,length(sel_set)-1,sel_variables,1)
    }else
    {
      sel_set <- which(as.matrix(coef(fit,s=0.07))[,1]!=0) #cv.glmnet select no variables, set a specific lambda
      if (length(sel_set)>1) # at least one variable selected
      {
        sel_variables=paste((sel_set-1)[2:length(sel_set)],collapse="," )
        Xsel= Xall[,(sel_set-1)[2:length(sel_set)]]
        fit1=lm(Y~Xsel)
        #summary(fit1)
        rsquared=summary(fit1)$r.squared
        Yhat=fit1$fitted.values
        res1[count,]=c(Yhat,rsquared,length(sel_set)-1,sel_variables,0)
      }
    }
    
    
    sel_set <- which(as.matrix(coef(fit,s=lambdas$lambda.1se))[,1]!=0)
    if (length(sel_set)>1) # at least one variable selected
    {
      sel_variables=paste((sel_set-1)[2:length(sel_set)],collapse="," )
      Xsel= Xall[,(sel_set-1)[2:length(sel_set)]]
      fit1=lm(Y~Xsel)
      #summary(fit1)
      rsquared=summary(fit1)$r.squared
      Yhat=fit1$fitted.values
      res2[count,]=c(Yhat,rsquared,length(sel_set)-1,sel_variables,1)
    }else
    {
      sel_set <- which(as.matrix(coef(fit,s=0.1))[,1]!=0) #cv.glmnet select no variables, set a specific lambda
      if (length(sel_set)>1) # at least one variable selected
      {
        sel_variables=paste((sel_set-1)[2:length(sel_set)],collapse="," )
        Xsel= Xall[,(sel_set-1)[2:length(sel_set)]]
        fit1=lm(Y~Xsel)
        #summary(fit1)
        rsquared=summary(fit1)$r.squared
        Yhat=fit1$fitted.values
        res2[count,]=c(Yhat,rsquared,length(sel_set)-1,sel_variables,0)
      }
    }
    if (count%%50==0) print(count)
    count=count+1
  #}
  result=list(resmin=res1,res1se=res2)
  return(result)
}


output1=paste0("./predictmrna/predictedmrnamin_all1.txt")
output2=paste0("./predictmrna/predictedmrna1se_all1.txt")

mpi.bcast.Robj2slave(predictrmrna)
mpi.bcast.Robj2slave(X)
mpi.bcast.Robj2slave(mrna)
mpi.bcast.Robj2slave(cv.glmnets)
mpi.bcast.Robj2slave(predictrmrna_1row)
mpi.bcast.cmd(library(glmnet))
mpi.bcast.cmd(require(methods))

res1 <- NULL
res2 <- NULL
nrun <- ceiling(nrow(mrna)/1000)
for (j in 1:nrun){
  cat(j,"..")
  if (j < nrun) cseq <- ((j-1)*1000+1):(j*1000)  else  cseq <- ((j-1)*1000+1):nrow(mrna)
  z=mrna[cseq,]
  res=mpi.parRapply(X=z,FUN=predictrmrna_1row,Xall=X,job.num=njob)
  tmp=matrix(unlist(res),byrow=TRUE,ncol=ncol(mrna)+4) #for each line in mrna, return 2 lines resmin and res1se
  res1=rbind(res1,tmp[seq(1,nrow(tmp),2),]) #odd row resmin
  res2=rbind(res2,tmp[seq(2,nrow(tmp),2),]) #even row res1se
  rownames(res1)[cseq]=rownames(res2)[cseq]=names(res)
  colnames(res1)=colnames(res2)=c(colnames(mrna),"rsquared","numvariables","selvariables","cvwork")
  write.table(res1,file=output1,col.names=T,row.names=T,sep="\t",quote=F)
  write.table(res2,file=output2,col.names=T,row.names=T,sep="\t",quote=F)
}
Sys.time()
mpi.close.Rslaves()
#mpi.close.Rslaves(dellog = FALSE)
mpi.quit()

#test=mpi.parRapply(X=z,FUN=predictrmrna_1row,Xall=X,job.num=njob)

#rownames(test1)=rownames(mrna)[1]
#colnames(test1)=c(colnames(mrna),"rsquared","numvariables","selvariables","cvwork")

cmd=paste0("salloc -t 1-10 -n ",njob," mpirun -n 1 ","./mpi_predictmrna.R ",njob)
print(cmd)

#salloc -t 1-1 -n 30 mpirun -n 1 R --interactive
