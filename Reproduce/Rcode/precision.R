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

# simulation for precision determination
########################################

ndat.k <- c(2:8)
nsims<-10000 #number of simulations

set.seed(1)
n<-length(ndat.k)
CV.est<-matrix(NA,nrow=n,ncol=length(unique(dPCRHIV$exp)))
CV.for.plot<-matrix(NA,nsims,n)

for(k in 1:n){
ndat<-ndat.k[k] #number of replicates
CV.hat<-matrix(NA,nsims,length(unique(dPCRHIV$exp)))
for(g in 1:nsims){
if(g%%200==0){print(g)}
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),1:2]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
for(i in 1:nrow(sim.dat)){
	CV.hat[g,i]<-sd(sim.res[((i-1)*ndat+1):(i*ndat),2])/mean(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
}
CV.for.plot[,k]<-CV.hat[,1]
CV.est[k,]<-apply(CV.hat,2,sd)/apply(CV.hat,2,mean)
}

# make a plot

CV.df <- data.frame(CV=array(CV.for.plot),replicates=c(rep(c(2:8),each=10000)))

ggplot(CV.df, aes(y=CV, x=factor(replicates))) +
  geom_boxplot(alpha = 0.2) +
  ylim(0, 0.13) + coord_flip() +
  theme_minimal() +
  ylab(expression(R^2)) + xlab("") +
  theme(legend.position="none")

# scaled plot for improved interpretability

CV.for.plot.perc <- t(t(CV.for.plot)/(sqrt(sim.dat[1,1])/sim.dat[1,2]))*100-100
CV.df.perc <- data.frame(CV=array(CV.for.plot.perc),replicates=c(rep(c(2:8),each=10000)))

ggplot(CV.df.perc, aes(y=CV, x=factor(replicates))) +
  geom_boxplot(alpha = 0.2) +
  ylim(-100, 350) + coord_flip() +
  theme_minimal() + geom_hline(yintercept=0,alpha=0.2,linetype=2) +
  ylab("Deviation from true CV (%)") + xlab("Replicates") +
  theme(legend.position="none")