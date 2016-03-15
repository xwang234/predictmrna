#!/usr/bin/env Rscript

setwd("/fh/fast/dai_j/CancerGenomics/Ovarian")
load("./mrna_copynumber_methylation_mutation.RData")


library(gdata)
clinic=read.xls("/fh/fast/dai_j/CancerGenomics/Ovarian/TCGAdata/TCGA-OV-Clinical-Table_S1.2.xlsx")

clinic <- clinic[order(clinic[,1]),]

mrnacluster <- mrnacluster[order(mrnacluster[,1]),] 

names(clinic)[1] <- "ID"

names(mrnacluster)[1] <- "ID"

clinic <- merge(clinic,mrnacluster,by="ID",x.all=T,y.all=F)

load("mrna.RData")
ID <- as.character(names(mrna))

list1 <- as.character(clinic[clinic[,7]=="PROGRESSIVE DISEASE"| clinic[,7]=="STABLE DISEASE",1])

list2 <- as.character(clinic[clinic[,13]=="Resistant" & clinic[,7]!="PROGRESSIVE DISEASE" & clinic[,7]!="STABLE DISEASE" ,1])

list3 <- as.character(clinic[clinic[,13]=="Sensitive" & !is.na(clinic[,12])& clinic[,12]<12& clinic[,7]!="PROGRESSIVE DISEASE" & clinic[,7]!="STABLE DISEASE",1])

list4 <- as.character(clinic[clinic[,13]=="Sensitive" & !is.na(clinic[,12])& clinic[,12]<=24 & clinic[,12]>=12& clinic[,7]!="PROGRESSIVE DISEASE" & clinic[,7]!="STABLE DISEASE",1])

list5 <-as.character(clinic[clinic[,13]=="Sensitive" & !is.na(clinic[,12])& clinic[,12]>24 & clinic[,7]!="PROGRESSIVE DISEASE" & clinic[,7]!="STABLE DISEASE",1,1])

group <- rep(NA,489)

for (i in 1:489) {
  if (sum(ID[i]==list1)==1) group[i] <- 0
  if (sum(ID[i]==list2)==1) group[i] <- 1
  if (sum(ID[i]==list3)==1) group[i] <- 2
  if (sum(ID[i]==list4)==1) group[i] <- 2
  if (sum(ID[i]==list5)==1) group[i] <- 2
}

for (i in 1:489) mrna[,i] <- as.numeric(mrna[,i])

mdat <- t(mrna)
for (i in 1:11864) mdat[,i] <- as.numeric(mdat[,i])
for (i in 1:11864) mdat[,i] <- mdat[,i]/sqrt(var(mdat[,i])) 



library(glmnet)
 x=mdat[!is.na(group),]
 y=group[!is.na(group)]
 fit <- glmnet(x,y,family="multinomial")
 cvfit <- cv.glmnet(x,y,family="multinomial")
 
 fit <- glmnet(x,y)
 cvfit <- cv.glmnet(x,y,nfolds=length(y))
 
test=coef(cvfit,s=cvfit$lambda.min)
test1=as.matrix(test[[1]])
test2=test1[test1[,1]>0,1]
test1=as.matrix(test[[2]])
test2=test1[test1[,1]>0,1]
test1=as.matrix(test[[3]])
test2=test1[test1[,1]>0,1]

glmnetclassify=function(data,group)
{
  #remove rows contain NA
  tmp=rowSums(data)
  idx=is.na(tmp)
  data=data[!idx,]
  mdat <- t(data)
  for (i in 1:nrow(data)) mdat[,i] <- as.numeric(mdat[,i])
  for (i in 1:nrow(data)) mdat[,i] <- mdat[,i]/sqrt(var(mdat[,i]))
  x=mdat[!is.na(group),]
  y=group[!is.na(group)]
  y=as.factor(y)
  fit <- glmnet(x,y,family="multinomial")
  cvfit <- cv.glmnet(x,y,family="multinomial",nfolds = length(y))
  pfit = as.integer(predict(fit,x,s=cvfit$lambda.min,type="class"))
  res=list(cvfit=cvfit,fit=fit,pfit=pfit,error=sum(pfit!=y)/length(y))
  return(res)
}
data=mrna
rescmrna=glmnetclassify(data,group)
predictmrna=predictmrna1se
data=predictmrna[,1:489]
rescpredictmrna1se=glmnetclassify(data,group)
predictmrna=predictmrnamin
data=predictmrna[,1:489]
rescpredictmrnamin=glmnetclassify(data,group)

data=t(X) #X from predictmrna.RData
resX=glmnetclassify(data,group)
data=rbind(mrna,t(X))
resmrnaX=glmnetclassify(data,group)

glmnetclassify2=function(data,group)
{
  #remove rows contain NA
  tmp=rowSums(data)
  idx=is.na(tmp)
  data=data[!idx,]
  idx=group==0 #remove group0
  group[idx]=NA
  mdat <- t(data)
  for (i in 1:nrow(data)) mdat[,i] <- as.numeric(mdat[,i])
  for (i in 1:nrow(data)) mdat[,i] <- mdat[,i]/sqrt(var(mdat[,i]))
  x=mdat[!is.na(group),]
  y=group[!is.na(group)]
  y=as.factor(y)
  fit <- glmnet(x,y,family="binomial")
  cvfit <- cv.glmnet(x,y,family="binomial",nfolds = length(y))
  pfit = as.integer(predict(fit,x,s=cvfit$lambda.min,type="class"))
  res=list(cvfit=cvfit,fit=fit,pfit=pfit,error=sum(pfit!=y)/length(y))
  return(res)
}
data=mrna
resc2mrna=glmnetclassify2(data,group)
predictmrna=predictmrna1se
data=predictmrna[,1:489]
resc2predictmrna1se=glmnetclassify2(data,group)
predictmrna=predictmrnamin
data=predictmrna[,1:489]
resc2predictmrnamin=glmnetclassify2(data,group)


# x=mdat[!is.na(group),]
# y=group[!is.na(group)]
# #fit <- glmnet(mdat[!is.na(group),],group[!is.na(group)])
# fit <- glmnet(x,y)
# cvfit <- cv.glmnet(mdat[!is.na(group),],group[!is.na(group)],nfolds=length(group[!is.na(group)]))
# 
# cvfit$lambda.1se
# cvfit$lambda.min
# 
# sel_set <- which(as.matrix(coef(fit,s=cvfit$lambda.1se))[,1]!=0)
# 
# pfit = predict(fit,mdat[!is.na(group),],s=cvfit$lambda.1se,type="response")
# pfit = predict(fit,mdat[!is.na(group),],s=cvfit$lambda.min,type="response")
# 
# xx <- cbind(pfit,group[!is.na(group)])
# #xx[,2] <- ifelse(xx[,2]>=1,1,0)
# xx[,2] <- ifelse(xx[,2]>=2,1,0)
# xx <- xx[order(xx[,1]),]
# ss <- xx
# for (i in 1:nrow(xx)) {
#   ss[i,1] <- mean(xx[xx[,2]==0,1]>=xx[i,1])
#   ss[i,2] <- mean(xx[xx[,2]==1,1]>=xx[i,1])
# }
# 
# sens=ss[,2]
# spec=ss[,1]
# idx=order(spec)
# spec=spec[idx]
# sens=sens[idx]
# height = (sens[-1]+sens[-length(sens)])/2
# width = diff(spec)
# sum(height*width)
# 
# plot(ss[,1],ss[,2],xlab="False positive rate",ylab="True positive rate",type="l",col=3,lwd=4)
# plot(xx)
col_rsquared=490 #column of r-squared
col_cvwork=493 #if cv.glmnet selected variables
col_numvar=491 #number of variables
col_variable=492 #selected variables
formcv=function(predictmrna)
{
  if (ncol(predictmrna)>489)
  {
    print(sum(is.na(predictmrna[,col_rsquared])))
    
    idxNA=is.na(predictmrna[,col_rsquared])
    predictmrna1=predictmrna[!idxNA,]
    #predictdmrna1[idxNA,1:ncol(mrna)]=mrna[idxNA,]
    hist(predictmrna[,col_rsquared],xlab="Rsquared",main="")
    hist(predictmrna[,col_numvar],xlab="Num of features",main="")
    mean(predictmrna[,col_rsquared],na.rm=T)
  }else
  {
    predictmrna1=predictmrna
  }
  
  #[1] 0.411332
  
  #for (i in 1:489) predictmrna1[,i] <- as.numeric(predictmrna1[,i])
  
  mdat1 <- t(predictmrna1[,1:489])
  for (i in 1:ncol(mdat1)) mdat1[,i] <- as.numeric(mdat1[,i])
  for (i in 1:ncol(mdat1)) mdat1[,i] <- mdat1[,i]/sqrt(var(mdat1[,i])) 
  
  library(glmnet)
  x1=mdat1[!is.na(group),]
  y=group[!is.na(group)]
  #fit1 <- glmnet(mdat1[!is.na(group),],group[!is.na(group)])
  fit1 <- glmnet(x1,y)
  cvfit1 <- cv.glmnet(mdat1[!is.na(group),],group[!is.na(group)],nfolds=length(group[!is.na(group)]))
  print(cvfit1$lambda.1se)
  print(cvfit1$lambda.min)
  sel_set1se <- which(as.matrix(coef(fit1,s=cvfit1$lambda.1se))[,1]!=0)
  sel_setmin = which(as.matrix(coef(fit1,s=cvfit1$lambda.min))[,1]!=0)
  pfit1se=pfitmin=rep(NA,length(y))
  pfitmin = predict(fit1,mdat1[!is.na(group),],s=cvfit1$lambda.min,type="response")
  pfit1se = predict(fit1,mdat1[!is.na(group),],s=cvfit1$lambda.1se,type="response")
  result=list(cv=cvfit1,fit=fit1,sel_setmin=sel_setmin,sel_set1se=sel_set1se,pfitmin=pfitmin,pfit1se=pfit1se)
}

plotroc=function(pfit)
{
  xx <- cbind(pfit,group[!is.na(group)])
  #xx[,2] <- ifelse(xx[,2]>=1,1,0)
  xx[,2] <- ifelse(xx[,2]>=2,1,0)
  xx <- xx[order(xx[,1]),]
  ss <- xx
  for (i in 1:nrow(xx)) {
    ss[i,1] <- mean(xx[xx[,2]==0,1]>=xx[i,1])
    ss[i,2] <- mean(xx[xx[,2]==1,1]>=xx[i,1])
  }
  
  sens=ss[,2]
  spec=ss[,1]
  idx=order(spec)
  spec=spec[idx]
  sens=sens[idx]
  height = (sens[-1]+sens[-length(sens)])/2
  width = diff(spec)
  print(sum(height*width))
  
  plot(ss[,1],ss[,2],xlab="False positive rate",ylab="True positive rate",type="l",col=3,lwd=4)
  
}


resmrna=formcv(mrna)
plotroc(resmrna$pfitmin)
plotroc(resmrna$pfit1se)


predictmrnamin=NULL
for (i in 1:30)
{
  tmp=read.table(file=paste0("./predictmrna/predictedmrnamin_",i,".txt"),header=T,sep="\t")
  predictmrnamin=rbind(predictmrnamin,tmp)
}

predictmrna1se=NULL
for (i in 1:30)
{
  tmp=read.table(file=paste0("./predictmrna/predictedmrna1se_",i,".txt"),header=T,sep="\t")
  predictmrna1se=rbind(predictmrna1se,tmp)
}

predictmrnamin=read.table(file="./predictmrna/predictedmrnamin_all2.txt",header=T,sep="\t")
predictmrna1se=read.table(file="./predictmrna/predictedmrna1se_all2.txt",header=T,sep="\t")
predictmrnamin=read.table(file="./predictmrna/predictedmrnamin_all4.txt",header=T,sep="\t")
predictmrna1se=read.table(file="./predictmrna/predictedmrna1se_all4.txt",header=T,sep="\t")

predictmrnamin[,col_variable]=as.character(predictmrnamin[,col_variable])
predictmrna1se[,col_variable]=as.character(predictmrna1se[,col_variable])

hist(predictmrnamin$numvariables)
hist(predictmrna1se$numvariables)

hist(predictmrnamin$rsquared)
hist(predictmrna1se$rsquared)

respredictmrnamin=formcv(predictmrnamin)
plotroc(respredictmrnamin$pfitmin)
plotroc(respredictmrnamin$pfit1se)

respredictmrna1se=formcv(predictmrna1se)
plotroc(respredictmrna1se$pfitmin)
plotroc(respredictmrna1se$pfit1se)

#combine 1se and min results
predictmrna=predictmrna1se
idx_cvworkNA=is.na(predictmrna[,col_cvwork])
predictmrna[idx_cvworkNA,]=predictmrnamin[idx_cvworkNA,]
hist(predictmrna$numvariables,main="",xlab="Number of features")
hist(predictmrna$rsquared,main="",xlab="R-squared")

resall1=formcv(predictmrna)
plotroc(resall$pfitmin)


maxnumvar=100
idxNA=!is.na(predictmrna$numvariables) & predictmrna$numvariables>=maxnumvar
predictmrna[idxNA,]=rep(NA,ncol(predictmrna))
resnum1=formcv(predictmrna)
plotroc(resnum1$pfitmin)

min(predictmrna$numvariables,na.rm=T)
min(predictmrna$rsquared,na.rm=T)
minrsquared=0.1
idxNA=!is.na(predictmrna$rsquared) & predictmrna$rsquared<=minrsquared
predictmrna[idxNA,]=rep(NA,ncol(predictmrna))
resrsquared=formcv(predictmrna)
plotroc(resrsquared$pfitmin)

resrsquarednum=formcv(predictmrna) #consider both variable num and rsquared
plotroc(resrsquarednum$pfitmin)

#save(resmrna,resall,resrsquared,resrsquarednum,resnum,file="analyze_predictmrna0229.RData")
test=names(resmrna$sel_set1se) %in% names(resmrna$sel_setmin)
sum(test)
#[1] 27

test=names(resmrna$sel_set1se) %in% names(resall$sel_setmin)
sum(test)
#[1] 9

features=colnames(X)

countseltypes=function(features,predictmrna)
{
  col_variable=492
  numtypes=c(10391,9329,5876)
  res=apply(predictmrna,1,function(x1) {
    result=rep(0,3)
    if (!is.na(x1[[col_variable]]))
    {
      tmp=as.numeric(unlist(strsplit(as.character(x1[[col_variable]]),",")))
      tmp1=sum(tmp<=numtypes[1])
      tmp2=sum(numtypes[1]<tmp & tmp<=numtypes[1]+numtypes[2])
      tmp3=sum(tmp>numtypes[1]+numtypes[2])
      if (tmp1>0) result[1]=tmp1
      if (tmp2>0) result[2]=tmp2
      if (tmp3>0) result[3]=tmp3
    }
    return(result)
  })
  res=data.frame(t(res))
  colnames(res)=c("CP","ME","MU")
  rownames(res)=rownames(predictmrna)
  return(res)
}

predictmrna=predictmrna1se
predictmrna=predictmrnamin
predictmrna=predictmrna[predictmrna$cvwork==1,]
count1se=countseltypes(features,predictmrna)
hist(count1se[,1],main="Copy number")
hist(count1se[,2],main="Methylation")
hist(count1se[,3],main="Mutation")
hist(rowSums(count1se),main="ALL")
count1se_all=rowSums(count1se)
sum(count1se_all==1)
sum(predictmrna$rsquared<0.1,na.rm=T)

cp_prop=apply(count1se,1,function(x){x[1]/sum(x)})
me_prop=apply(count1se,1,function(x){x[2]/sum(x)})
mu_prop=apply(count1se,1,function(x){x[3]/sum(x)})
hist(cp_prop,main="Copy number")
hist(me_prop,main="Methylation")
hist(mu_prop,main="Mutation")
data=data.frame(matrix(NA,nrow=3,ncol=nrow(cor_mrna_copynumber)))
colnames(data)=rownames(cor_mrna_copynumber)
data[1,]=cp_prop
data[2,]=me_prop
data[3,]=mu_prop
colors=c("red","blue","green")
barplot(as.matrix(data[,1:100]), main="", ylab = "Proportion", cex.lab = 1, cex.main = 1, beside=TRUE, col=colors)
legend("topleft", c("CP","ME","MU"), cex=1.3, bty="n", fill=colors)





sum(is.na(predictmrna[,col_rsquared]))
# [1] 170
idxNA=is.na(predictmrna[,col_rsquared])
predictmrna1=predictmrna[!idxNA,]
#predictdmrna1[idxNA,1:ncol(mrna)]=mrna[idxNA,]
hist(predictmrna[,col_rsquared],xlab="Rsquared",main="")
hist(predictmrna[,col_numvar],xlab="Num of features",main="")
mean(predictmrna[,col_rsquared],na.rm=T)
#[1] 0.411332

for (i in 1:489) predictmrna1[,i] <- as.numeric(predictmrna1[,i])

mdat1 <- t(predictmrna1[,1:489])
for (i in 1:ncol(mdat1)) mdat1[,i] <- as.numeric(mdat1[,i])
for (i in 1:ncol(mdat1)) mdat1[,i] <- mdat1[,i]/sqrt(var(mdat1[,i])) 

library(glmnet)
x1=mdat1[!is.na(group),]
y=group[!is.na(group)]
#fit1 <- glmnet(mdat1[!is.na(group),],group[!is.na(group)])
fit1 <- glmnet(x1,y)

cvfit1 <- cv.glmnet(x1,y)
cvfit1 <- cv.glmnet(mdat1[!is.na(group),],group[!is.na(group)],nfolds=length(group[!is.na(group)]))
#cvfit1 <- cv.glmnet(mdat[!is.na(group),],group[!is.na(group)],nfolds=length(group[!is.na(group)]))
cvfit1$lambda.1se
cvfit1$lambda.min

fit2=glmnet(x1,y,family="multinomial")
cvfit2=cv.glmnet(x1,y,family="multinomial",nfolds=length(group[!is.na(group)]))
cvfit2$lambda.1se
cvfit2$lambda.min

sel_set1 <- which(as.matrix(coef(fit2,s=cvfit2$lambda.min))[,1]!=0)

sel_set1 <- which(as.matrix(coef(fit1,s=cvfit1$lambda.min))[,1]!=0)
sel_set1 <- which(as.matrix(coef(fit1,s=cvfit1$lambda.1se))[,1]!=0)

pfit = predict(fit1,mdat1[!is.na(group),],s=cvfit1$lambda.min,type="response")

sel_set2 <- which(as.matrix(coef(fit1keep,s=cvkeep$lambda.min))[,1]!=0)



#use lambda.min
predictdmrnamin=NULL
for (i in 1:30)
{
  tmp=read.table(file=paste0("./predictmrna/predictedmrnamin_",i,".txt"),header=T,sep="\t")
  #tmp=read.table(file=paste0("./predictmrna/predictedmrna1se_",i,".txt"),header=T,sep="\t")
  predictdmrnamin=rbind(predictdmrnamin,tmp)
}
sum(is.na(predictdmrnamin[,col_rsquared]))
# [1] 15
idxNA=is.na(predictdmrnamin[,col_rsquared])
predictdmrnamin1=predictdmrnamin[!idxNA,]
#predictdmrna1[idxNA,1:ncol(mrna)]=mrna[idxNA,]
hist(predictdmrnamin[,col_rsquared],xlab="Rsquared",main="Lambda.min")
mean(predictdmrnamin[,col_rsquared],na.rm=T)
mdat2 <- t(predictdmrnamin1[,1:489])
for (i in 1:ncol(mdat2)) mdat2[,i] <- as.numeric(mdat2[,i])
for (i in 1:ncol(mdat2)) mdat2[,i] <- mdat2[,i]/sqrt(var(mdat2[,i])) 


x2=mdat2[!is.na(group),]
y=group[!is.na(group)]
#fit1 <- glmnet(mdat1[!is.na(group),],group[!is.na(group)])
fit2 <- glmnet(x2,y)
cvfit2 <- cv.glmnet(x2,y)
cvfit2 <- cv.glmnet(x2,y,nfolds=length(y))
#cvfit1 <- cv.glmnet(mdat[!is.na(group),],group[!is.na(group)],nfolds=length(group[!is.na(group)]))
cvfit2$lambda.1se
cvfit2$lambda.min
sel_set2 <- which(as.matrix(coef(fit2,s=cvfit2$lambda.min))[,1]!=0)
pfit = predict(fit2,mdat2[!is.na(group),],s=cvfit2$lambda.min,type="response")

#check predicted mrna
load("predictmrna.RData") #contains the X matrix
gene="ABCA12"
gene="A4GNT" # 1 variable selected in 1se
gene="ABCB4" # 1se has no variable selected
idx=which(rownames(mrna)==gene)
Ytest=t(as.matrix(mrna[idx,]))
# Y=Y[!idx]
fittest = glmnet(X, Ytest)
cvfittest=cv.glmnet(X,Ytest,nfolds=length(Ytest))
cvfittest$lambda.min
cvfittest$lambda.1se
plot(cvfittest)

cvfittest1=cv.glmnet(X,Ytest)
cvfittest1$lambda.min
cvfittest1$lambda.1se
lambdas =cv.glmnets(X,Ytest,nfolds=10,ntime=10)

sel_settest <- which(as.matrix(coef(fittest,s=lambdas$lambda.1se))[,1]!=0)
sel_settest <- which(as.matrix(coef(fittest,s=cvfittest1$lambda.1se))[,1]!=0)

cvfittest2=cv.glmnet(X,Ytest)
cvfittest2$lambda.min
cvfittest2$lambda.1se
lambdas =cv.glmnets(X,Ytest,nfolds=10,ntime=10)
