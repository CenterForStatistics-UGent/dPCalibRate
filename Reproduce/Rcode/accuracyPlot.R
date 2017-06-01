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





#accuracy plot for paper
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=3)
set.seed(1)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),1:2]<-cbind(rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])),sim.dat[i,2])
}

sim.res[,3]<-sim.res[,2]
dPCRHIV<-data.frame(sim.res)
dPCRHIV[1:ndat,1]<-dPCRHIV[1:ndat,1]*0.9 #what happens when we introduce a bias?
colnames(dPCRHIV)<-c("obs","exp")
a <- lm(sim.res[,3]~0+sim.res[,2])
#preds <- unique(dPCRHIV$exp)
preds <-unique(predict(a, newdata=data.frame(exp=unique(dPCRHIV$exp))))
groupmeans<-rowMeans(matrix(dPCRHIV$obs,ncol=5,byrow=T))

plot <- list()
for(i in 1:length(preds)){
dataSubs <- data.frame(type=c("pred",rep("obs",ndat)), data=c(preds[i],dPCRHIV$obs[c(((i-1)*ndat+1):(i*ndat))]),exp=unique(dPCRHIV$exp)[i])
if(max(dataSubs$data[-1])/dataSubs$data[1]>1/(min(dataSubs$data[-c(1,which(dataSubs$data==0))])/dataSubs$data[1])){
	cut <- (floor(max(dataSubs$data[-1])/dataSubs$data[1]*100)-100)/2/100
	level<-1+cut*2
	yint<- c()
	while(level>min(dataSubs$data[-1])/dataSubs$data[1]){
		if(level!=1){yint<-c(yint,level)}
		level<-level-cut
	}
} else {
	cut <- (ceiling(min(dataSubs$data[-1])/dataSubs$data[1]*100)-100)/2/100
	level<-1+cut*2
	yint<- c()
	while(level<max(dataSubs$data[-1])/dataSubs$data[1]){
		if(level!=1){yint<-c(yint,level)}
		level<-level-cut
	}
}

yint<-yint*unique(dataSubs$data[1])
if(length(yint)==0){
yint<-0
}
ylabs<-paste(round((yint/unique(dataSubs$data[1])-1)*100,1),"%",sep="")
if(length(which(ylabs=="0%"))>0){
yint<-yint[-which(ylabs=="0%")]
ylabs<-ylabs[-which(ylabs=="0%")]
}

if(i==7){
plot[[i]] <- ggplot(dataSubs, aes(x=exp,y=data,shape=type)) +
  stat_smooth(method="lm",se=FALSE,size=0.2) +  
  theme_minimal() + geom_hline(yintercept=preds[i]) +
  geom_hline(yintercept=groupmeans[i],linetype=2) + geom_point() +
  geom_hline(yintercept=yint,linetype=3,alpha=0.25) + annotate("text",label=ylabs,y=yint,x=unique(dataSubs$exp),alpha=0.75,size=3) +
  scale_colour_brewer(palette = "Set1") +
  xlab("")+ylab("observed concentration") + theme(legend.position="none", axis.ticks.y=element_blank(),axis.text.y=element_blank()) + scale_x_discrete(limits=c(unique(dPCRHIV$exp)[i]))	
} else if(i==4){
	plot[[i]] <- ggplot(dataSubs, aes(x=exp,y=data,shape=type)) +
  stat_smooth(method="lm",se=FALSE,size=0.2) +  
  theme_minimal() + geom_hline(yintercept=preds[i]) +
  geom_hline(yintercept=groupmeans[i],linetype=2) + geom_point() +
  geom_hline(yintercept=yint,linetype=3,alpha=0.25) + annotate("text",label=ylabs,y=yint,x=unique(dataSubs$exp),alpha=0.75,size=3) +
  scale_colour_brewer(palette = "Set1") +
  xlab("expected concentration")+ylab("") + theme(legend.position="none", axis.ticks.y=element_blank(),axis.text.y=element_blank()) + scale_x_discrete(limits=c(unique(dPCRHIV$exp)[i]))
} else {
plot[[i]] <- ggplot(dataSubs, aes(x=exp,y=data,shape=type)) +
  stat_smooth(method="lm",se=FALSE,size=0.2) +  
  theme_minimal() + geom_hline(yintercept=preds[i]) +
  geom_hline(yintercept=groupmeans[i],linetype=2) + geom_point() +
  geom_hline(yintercept=yint,linetype=3,alpha=0.25) + annotate("text",label=ylabs,y=yint,x=unique(dataSubs$exp),alpha=0.75,size=3) +
  scale_colour_brewer(palette = "Set1") +
  xlab("")+ylab("") + theme(legend.position="none", axis.ticks.y=element_blank(),axis.text.y=element_blank()) + scale_x_discrete(limits=c(unique(dPCRHIV$exp)[i]))
}
}
#reorder
plot2<-plot
for(i in 1:7){
	plot2[[i]]<-plot[[8-i]]
}
plot_grid(plotlist=plot2, align='h', labels=c('', '', '', '', '', '', ''),nrow=1)