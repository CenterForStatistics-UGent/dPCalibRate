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



# dPCR case study data preparation
##################################

dPCRHIV<-c(t(data[17:24,3:7]))
dPCRHIV<-data.frame(cbind(dPCRHIV,sort(rep(data$V2[1:8],5),decreasing=T)))
dPCRHIV<-dPCRHIV[-c(36:40),]
colnames(dPCRHIV)<-c("obs","exp")
CV.CI<-list()
for(i in 1:7){
	CV.CI[[i]]<-ci.cv(data=dPCRHIV[((i-1)*5+1):(i*5),1])
}
var.calc<-NA
for(i in 1:7){
	var.calc[i]<-var(dPCRHIV[((i-1)*5+1):(i*5),1])
}
sim.dat<-data.frame(var.calc,unique(dPCRHIV$exp))



# analyse case study data
#########################

# quadratic regression and lack of fit test
ndat<-5
weights<-1/rep(var.calc,each=ndat)
fit<-lm(sim.res[,2]~sim.res[,1],weights=weights)
summary(lm(obs~exp,data=dPCRHIV))
summary(lm(obs~exp,data=dPCRHIV,weights=weights))
anova(lm(obs~(exp),data=dPCRHIV,weights=weights),lm(obs~factor(exp),data=dPCRHIV,weights=weights))
summary(lm(obs~exp+I(exp^2),data=dPCRHIV))
summary(lm(obs~exp+I(exp^2),data=dPCRHIV,weights=weights))

# frequency test
N <- length(unique(dPCRHIV$exp))
vec <- lm(sim.res[,2]~sim.res[,1],weights=weights)$residuals
1-pgamma(sum((colSums(matrix(vec>0,ncol=7))-ndat*0.5)^2/(ndat*0.5)),N/2)

# runs test
library(randtests)
res <-summary(lm(obs~exp,data=dPCRHIV,weights=weights))$residuals
runs.test(res,threshold=0,pvalue="exact")
res <-summary(lm(obs~exp,data=dPCRHIV))$residuals
runs.test(res,threshold=0,pvalue="exact")



# create case study linearity and accuracy plot
###############################################

a <- lm(obs~exp,data=dPCRHIV,weights=weights) #normal linear model
#preds <- unique(dPCRHIV$exp)
preds <-predict(a, newdata=data.frame(exp=unique(dPCRHIV$exp)))
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
plot_grid(plotlist=plot, align='h', labels=c('', '', '', '', '', '', '','B', '', '', '', '', '', ''),nrow=1)


plot2<-plot

dPCRHIV$exp1<-dPCRHIV$exp
a <- lm(exp1~0+exp,data=dPCRHIV,weights=weights) #normal linear model
#preds <- unique(dPCRHIV$exp)
preds <-predict(a, newdata=data.frame(exp=unique(dPCRHIV$exp)))
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
plot3<-plot2
for(i in 1:7){
	plot3[[i]]<-plot2[[8-i]]
}
plot4<-plot
for(i in 1:7){
	plot4[[i]]<-plot[[8-i]]
}
plot_grid(plotlist=c(plot3,plot4), align='h', labels=c('A', '', '', '', '', '', '','B', '', '', '', '', '', ''),nrow=2)