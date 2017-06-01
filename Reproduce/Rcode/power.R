# load libraries and data
#########################
library(MASS)
library(ggplot2)
library(cowplot)
library(sandwich)
library(lmtest)
library(reshape2)

data<-read.csv("Jones2016_BDQ.csv",header=F)
data$means<-rowMeans(data[,3:7])


# dPCR experiment 1 data preparation
####################################
dPCRHIV<-c(t(data[1:8,3:7])) #only data for first dPCR experiment
dPCRHIV<-data.frame(cbind(dPCRHIV,sort(rep(data$V2[1:8],5),decreasing=T)))
dPCRHIV<-dPCRHIV[-c(36:40),] #remove NTCs
colnames(dPCRHIV)<-c("obs","exp")

# calculate the variance at different concentrations
var.calc<-NA
for(i in 1:7){
	var.calc[i]<-var(dPCRHIV[((i-1)*5+1):(i*5),1])
}
sim.dat<-data.frame(var.calc,unique(dPCRHIV$exp))






# Calculate power for different tests
#####################################

# Quadratic regression

bias.vec<-c(0,0.10,0.20,0.20,0.20)
bias2.vec<-c(0,0,0,0.10,0.20)
ndat.k <- 2:8
set.seed(1)
power <- matrix(NA,length(ndat.k),length(bias.vec))
power.noweights <- matrix(NA,length(ndat.k),length(bias.vec))
power.hccm <- matrix(NA,length(ndat.k),length(bias.vec))

n<-length(ndat.k)
for(l in 1:length(bias.vec)){
for(k in 1:length(ndat.k)){
sim.dat
ndat<-ndat.k[k] #number of replicates
nsims<-10000 #number of simulations
bias<-bias.vec[l] #introduce a bias
bias2<-bias2.vec[l]
rsq<-array(NA,nsims)
rsq.noweights<-array(NA,nsims)
rsq.hccm<-array(NA,nsims)

for(g in 1:nsims){
if(g%%200==0){print(g)}
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),1:2]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
sim.res[(ndat+1):(2*(ndat)),2]<-sim.res[(ndat+1):(2*(ndat)),2]*(1-bias2)
var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
#WLS
weights<-1/rep(var.calc,each=ndat)
rsq[g]<-summary(lm(sim.res[,2]~sim.res[,1]+I(sim.res[,1]^2),weights=weights))$coefficients[3,4]

#OLS
rsq.noweights[g]<-summary(lm(sim.res[,2]~sim.res[,1]+I(sim.res[,1]^2)))$coefficients[3,4]

#WLS HC3
fit<-lm(sim.res[,2]~sim.res[,1]+I(sim.res[,1]^2),weights=weights)
rsq.hccm[g]<-coeftest(fit,vcov=vcovHC(fit,type="HC3"))[3,4]

}
print(c(l,k))
power[k,l]<-mean(rsq<0.05)
power.noweights[k,l]<-mean(rsq.noweights<0.05,na.rm=T)
power.hccm[k,l]<-mean(rsq.hccm<0.05,na.rm=T)
}
}
save.image("power.RData")


typeI<-data.frame(OLS=power.noweights[,1],WLS=power[,1],WLS.HCCM=power.hccm[,1])
typeI<-melt(typeI)
typeI<-data.frame(typeI,replicates=c(rep(2:8,3)))
colnames(typeI)<-c("type","size","replicates")
ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05)






# Runs test

library(randtests)
bias.vec<-c(0,0.10,0.20,0.20,0.20)
bias2.vec<-c(0,0,0,0.10,0.20)
ndat.k <- 2:8
set.seed(1)
power <- matrix(NA,length(ndat.k),length(bias.vec))
power.noweights <- matrix(NA,length(ndat.k),length(bias.vec))
power.hccm <- matrix(NA,length(ndat.k),length(bias.vec))
n<-length(ndat.k)
for(l in 1:length(bias.vec)){
for(k in 1:length(ndat.k)){
print(k)
ndat<-ndat.k[k] #number of replicates
nsims<-10000 #number of simulations
bias<-bias.vec[l] #introduce a bias
bias2<-bias2.vec[l]
rsq<-array(NA,nsims)
rsq.noweights<-array(NA,nsims)
rsq.hccm<-array(NA,nsims)

for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),1:2]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
sim.res[(ndat+1):(2*(ndat)),2]<-sim.res[(ndat+1):(2*(ndat)),2]*(1-bias2)
var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
#WLS
weights<-1/rep(var.calc,each=ndat)
fit<-lm(sim.res[,2]~sim.res[,1],weights=weights)
rsq[g]<-runs.test(fit$residuals,threshold=0,pvalue="exact")$p.value

#OLS
fit<-lm(sim.res[,2]~sim.res[,1])
rsq.noweights[g]<-runs.test(fit$residuals,threshold=0,pvalue="exact")$p.value


#WLS HC3
fit<-lm(sim.res[,2]~sim.res[,1],weights=weights)
rsq.hccm[g]<-runs.test(fit$residuals,threshold=0,pvalue="exact")$p.value
}
print(c(l,k))
power[k,l]<-mean(rsq<0.05)
power.noweights[k,l]<-mean(rsq.noweights<0.05,na.rm=T)
power.hccm[k,l]<-mean(rsq.hccm<0.05,na.rm=T)
}
}
save.image("powerruns.RData")
typeI<-data.frame(OLS=power.noweights[,1],WLS=power[,1],WLS.HCCM=power.hccm[,1])
typeI<-melt(typeI)
typeI<-data.frame(typeI,replicates=c(rep(2:8,3)))
colnames(typeI)<-c("type","size","replicates")
ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05)+ylim(0,0.8)



# Accuracy

bias.vec<-c(0,0.05,0.10,0.15,0.20)
ndat.k <- 2:8
set.seed(1)
power <- matrix(NA,length(ndat.k),length(bias.vec))
power.noweights <- matrix(NA,length(ndat.k),length(bias.vec))
power.hccm <- matrix(NA,length(ndat.k),length(bias.vec))
n<-length(ndat.k)
for(l in 1:length(bias.vec)){
for(k in 1:length(ndat.k)){
print(k)
ndat<-ndat.k[k] #number of replicates
nsims<-10000 #number of simulations
bias<-bias.vec[l] #introduce a bias
rsq<-array(NA,nsims)
rsq.noweights<-array(NA,nsims)
rsq.hccm<-array(NA,nsims)

for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),1:2]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[,2]<-sim.res[,2]*(1-bias)

var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
weights<-1/rep(var.calc,each=ndat)

#WLS
rsq[g]<-(1-pt(abs((summary(lm(sim.res[,2]~sim.res[,1],weights=weights))$coefficients[2,1]-1)/summary(lm(sim.res[,2]~sim.res[,1],weights=weights))$coefficients[2,2]),nrow(sim.res)-2))*2

#OLS
rsq.noweights[g]<-(1-pt(abs((summary(lm(sim.res[,2]~sim.res[,1]))$coefficients[2,1]-1)/summary(lm(sim.res[,2]~sim.res[,1]))$coefficients[2,2]),nrow(sim.res)-2))*2

#WLS HC3
fit<-lm(sim.res[,2]~sim.res[,1],weights=weights)
fit.se<-sqrt(vcovHC(fit,type="HC3")[2,2])
rsq.hccm[g]<-(1-pt(abs((summary(lm(sim.res[,2]~sim.res[,1],weights=weights))$coefficients[2,1]-1)/fit.se),nrow(sim.res)-2))*2
}
power[k,l]<-mean(rsq<0.05)
power.noweights[k,l]<-mean(rsq.noweights<0.05)
power.hccm[k,l]<-mean(rsq.hccm<0.05)

}
}
save.image("poweraccuracy.RData")



typeI<-data.frame(OLS=power.noweights[,1],WLS=power[,1],WLS.HCCM=power.hccm[,1])
typeI<-melt(typeI)
typeI<-data.frame(typeI,replicates=c(rep(2:8,3)))
colnames(typeI)<-c("type","size","replicates")
#ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05)+ylim(0,0.65)




# Frequency within a block test

bias.vec<-c(0,0.10,0.20,0.20,0.20)
bias2.vec<-c(0,0,0,0.10,0.20)
ndat.k <- c(2:8)
set.seed(1)
power <- matrix(NA,length(ndat.k),length(bias.vec))
power.noweights <- matrix(NA,length(ndat.k),length(bias.vec))
power.hccm <- matrix(NA,length(ndat.k),length(bias.vec))
n<-length(ndat.k)
for(l in 1:length(bias.vec)){
for(k in 1:length(ndat.k)){
print(k)
ndat<-ndat.k[k] #number of replicates
nsims<-10000 #number of simulations
bias<-bias.vec[l] #introduce a bias
bias2<-bias2.vec[l]
rsq<-array(NA,nsims)
rsq.noweights<-array(NA,nsims)
rsq.hccm<-array(NA,nsims)

for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),1:2]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
sim.res[(ndat+1):(2*(ndat)),2]<-sim.res[(ndat+1):(2*(ndat)),2]*(1-bias2)
var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}

N <- length(unique(dPCRHIV$exp))

weights<-1/rep(var.calc,each=ndat)
fit<-lm(sim.res[,2]~sim.res[,1],weights=weights)

vec <- fit$residuals
rsq[g]<-1-pgamma(sum((colSums(matrix(vec>0,ncol=7))-ndat*0.5)^2/(ndat*0.5)),N/2)
rsq.hccm[g]<-rsq[g]

vec <- lm(sim.res[,2]~sim.res[,1])$residuals
rsq.noweights[g]<-1-pgamma(sum((colSums(matrix(vec>0,ncol=7))-ndat*0.5)^2/(ndat*0.5)),N/2)


}
power[k,l]<-mean(rsq<0.05)
power.noweights[k,l]<-mean(rsq.noweights<0.05,na.rm=T)
power.hccm[k,l]<-mean(rsq.hccm<0.05,na.rm=T)

}
}
typeI<-data.frame(OLS=power.noweights[,1],WLS=power[,1],WLS.HCCM=power.hccm[,1],WLS.HCCM.aMINQUE=power.hccmMIN[,1],WLS.HCCM.aMINQUE.HAT=power.hccmMINhat[,1],IWLS.HCCM=power.irls.hccm[,1])
typeI<-melt(typeI)
typeI<-data.frame(typeI,replicates=c(rep(2:8,6)))
colnames(typeI)<-c("type","size","replicates")
ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05)+ylim(0,1)

save.image("powerfreq.RData")




# F test for lack of fit

bias.vec<-c(0,0.10,0.20,0.20,0.20)
bias2.vec<-c(0,0,0,0.10,0.20)
ndat.k <- c(2:8)
set.seed(1)
power <- matrix(NA,length(ndat.k),length(bias.vec))
power.noweights <- matrix(NA,length(ndat.k),length(bias.vec))
n<-length(ndat.k)
for(l in 1:length(bias.vec)){
for(k in 1:length(ndat.k)){
sim.dat
ndat<-ndat.k[k] #number of replicates
nsims<-10000 #number of simulations
bias<-bias.vec[l] #introduce a bias
bias2<-bias2.vec[l]
rsq<-array(NA,nsims)
rsq.noweights<-array(NA,nsims)
for(g in 1:nsims){
if(g%%200==0){print(g)}
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),1:2]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
sim.res[(ndat+1):(2*(ndat)),2]<-sim.res[(ndat+1):(2*(ndat)),2]*(1-bias2)
var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
weights<-1/rep(var.calc,each=ndat)

m1<-(lm(sim.res[,2]~as.factor(sim.res[,1]),weights=weights))
m2<-(lm(sim.res[,2]~(sim.res[,1]),weights=weights))
SSPE<-sum(weights*(m1$resid-predict(m1))^2)
SSPE.df<-nrow(sim.res)-length(unique(dPCRHIV$exp))
SSE<-sum(weights*(m2$resid-predict(m2))^2)
SSE.df<-nrow(sim.res)-2
F.w<-((SSE-SSPE)/(SSE.df-SSPE.df))/(SSPE/SSPE.df)
rsq[g]<-anova(m1,m2)[["Pr(>F)"]][2]
m1<-(lm(sim.res[,2]~as.factor(sim.res[,1])))
m2<-(lm(sim.res[,2]~(sim.res[,1])))
SSPE<-sum((m1$resid-predict(m1))^2)
SSPE.df<-nrow(sim.res)-length(unique(dPCRHIV$exp))
SSE<-sum((m2$resid-predict(m2))^2)
SSE.df<-nrow(sim.res)-2
F.w<-((SSE-SSPE)/(SSE.df-SSPE.df))/(SSPE/SSPE.df)
rsq.noweights[g]<-anova(m1,m2)[["Pr(>F)"]][2]
}
print(c(l,k))
power[k,l]<-mean(rsq<0.05)
power.noweights[k,l]<-mean(rsq.noweights<0.05,na.rm=T)
}
}
save.image("lackoffit.RData")