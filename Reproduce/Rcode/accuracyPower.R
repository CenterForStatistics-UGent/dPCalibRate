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






# Calculate power for accuracy
##############################

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


