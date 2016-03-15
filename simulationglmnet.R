data=data.frame(matrix(0,ncol=489,nrow=10))
for (i in 1:10)
{
  set.seed(i)
  x=rnorm(489,0.1*i,1+0.1*i)
  x=x/sd(x)
  x=x-min(x)
  x=x/max(x)
  data[i,]=x
 }
rownames(data)=paste0("feature",1:10)
set.seed(1)
beta=runif(10,-1,1)
beta
set.seed(1000)
y1=t(data)%*%beta
y2=0*rnorm(489)
y=y1+y2
hist(y1,main="target")
hist(y2,main="noise")
colnames(data)=colnames(mrna)
x1=rbind(mrna,data)
x1=t(x1)
library(glmnet)
fit1=glmnet(x1,y)
cvfit1=cv.glmnet(x1,y)
sel_set1 <- which(as.matrix(coef(fit1,s=cvfit1$lambda.min))[,1]!=0)
sel_set2 <- which(as.matrix(coef(fit1,s=cvfit1$lambda.1se))[,1]!=0)

cor2data=function(data1,data2)
{
  res=data.frame(matrix(NA,nrow=nrow(data1),ncol=nrow(data2)))
  rownames(res)=rownames(data1)
  colnames(res)=rownames(data2)
  for (i in 1:nrow(data2))
  {
    res[,i]=apply(data1,1,corxy,y=data2[i,])
  }
  return(res)
}

Xsel= x1[,(sel_set1-1)[2:length(sel_set1)]]
data1=t(Xsel)
data2=t(y)
res=cor2data(data1,data2)
plot(res)
res1=cor2data(t(x1),data2)
hist(as.numeric(res1[,1]))
as.numeric(max(res[,1]))

#cor of selected data and missing features
Xsel= x1[,(sel_set1-1)[2:length(sel_set1)]]
data1=t(Xsel)
missingfeatures=paste0("feature",c(1,2,3,8,9))
idx=colnames(x1) %in% missingfeatures
data2=t(x1[,idx])
res=cor2data(data1,data2)
testfit=lm(data2[3,]~t(data1))

data1=data
data2=t(y)
res=cor2data(data1,data2)

#cor of selected other features with selected real features
Xsel= x1[,(sel_set2-1)[2:length(sel_set2)]]
data1=t(Xsel)
idx=grepl("feature",rownames(data1))
data2=data1[!idx,]
data1=data1[idx,]
res=cor2data(data1,data2)
#don't have big cor