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



# simulation study for R^2 adequacy, OLS
########################################

set.seed(1)
sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.00 #introduce a bias
rsq<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
rsq[g]<-summary(lm(sim.res[,2]~sim.res[,1]))$r.squared
}

plot(density(rsq),xlim=c(0.98,1),main="",xlab=expression(R^2))
sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.10 #introduce a bias
rsq.bias<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
rsq.bias[g]<-summary(lm(sim.res[,2]~sim.res[,1]))$r.squared
}
lines(density(rsq.bias),col="blue")


sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.20 #introduce a bias
rsq.bias2<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
rsq.bias2[g]<-summary(lm(sim.res[,2]~sim.res[,1]))$r.squared
}
lines(density(rsq.bias2),col="green")



sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.20 #introduce a bias
bias2<-0.10 #introduce a second bias
rsq.bias4<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	

}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
sim.res[(ndat+1):(2*ndat),2]<-sim.res[(ndat+1):(2*ndat),2]*(1-bias2)
rsq.bias4[g]<-summary(lm(sim.res[,2]~sim.res[,1]))$r.squared
}
lines(density(rsq.bias4),col="grey")




sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.20 #introduce a bias
bias2<-0.20 #introduce a second bias
rsq.bias5<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	

}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
sim.res[(ndat+1):(2*ndat),2]<-sim.res[(ndat+1):(2*ndat),2]*(1-bias2)
rsq.bias5[g]<-summary(lm(sim.res[,2]~sim.res[,1]))$r.squared
}
lines(density(rsq.bias5),col="red")



# combined plot

bias.df1 <- data.frame(rsq=c(rsq,rsq.bias,rsq.bias2,rsq.bias4,rsq.bias5),bias=c(rep("No bias",10000),rep("Single bias 10%",10000),rep("Single bias 20%",10000),rep("Double bias 20% + 10%",10000),rep("Double bias 20% + 20%",10000)))
bias.df1$bias<-factor(bias.df1$bias, levels = levels(bias.df1$bias)[c(3,4,5,1,2)])
library(ggplot2)
library(RColorBrewer)
ggplot(bias.df1, aes(rsq, fill = bias, colour = bias)) +
  geom_density(alpha = 0.2, aes(linetype=bias)) +
  xlim(0.995, 1) +
  theme_minimal() +
  xlab(expression(R^2)) +
  theme(legend.title=element_blank()) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")

ggplot(bias.df1, aes(y=rsq, x=factor(bias))) +
  geom_boxplot(alpha = 0.2) +
  ylim(0.993, 1) + coord_flip() +
  theme_minimal() +
  ylab(expression(R^2)) + xlab("") +
  theme(legend.position="none")







# simulation study for R^2 adequacy, WLS
########################################

set.seed(1)
sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.00 #introduce a bias
rsq<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
a <- lm(sim.res[,2]~sim.res[,1])
weights<-1/rep(var.calc,each=ndat)

rsq[g]<-summary(lm(sim.res[,2]~sim.res[,1],weights=weights))$r.squared
}

 plot(density(rsq),xlim=c(0.98,1),main="",xlab=expression(R^2))
sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.10 #introduce a bias
rsq.bias<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
a <- lm(sim.res[,2]~sim.res[,1])
weights<-1/rep(var.calc,each=ndat)

rsq.bias[g]<-summary(lm(sim.res[,2]~sim.res[,1],weights=weights))$r.squared
}
lines(density(rsq.bias),col="blue")


sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.20 #introduce a bias
rsq.bias2<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
a <- lm(sim.res[,2]~sim.res[,1])
weights<-1/rep(var.calc,each=ndat)

rsq.bias2[g]<-summary(lm(sim.res[,2]~sim.res[,1],weights=weights))$r.squared
}
lines(density(rsq.bias2),col="green")




sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.20 #introduce a bias
bias2<-0.10 #introduce a second bias
bias3<-0.0
rsq.bias4<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
sim.res[(ndat+1):(2*ndat),2]<-sim.res[(ndat+1):(2*ndat),2]*(1-bias2)
var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
a <- lm(sim.res[,2]~sim.res[,1])
weights<-1/rep(var.calc,each=ndat)

rsq.bias4[g]<-summary(lm(sim.res[,2]~sim.res[,1],weights=weights))$r.squared
}
lines(density(rsq.bias4),col="grey")




sim.dat
ndat<-5 #number of replicates
nsims<-10000 #number of simulations
bias<-0.20 #introduce a bias
bias2<-0.20 #introduce a second bias
rsq.bias5<-array(NA,nsims)
for(g in 1:nsims){
sim.res<-matrix(NA,nrow=nrow(sim.dat)*ndat,ncol=2)
for(i in 1:nrow(sim.dat)){
	sim.res[c(((i-1)*ndat+1):(i*ndat)),]<-cbind(sim.dat[i,2],rnorm(ndat,sim.dat[i,2],sqrt(sim.dat[i,1])))	
}
sim.res[1:ndat,2]<-sim.res[1:ndat,2]*(1-bias)
sim.res[(ndat+1):(2*ndat),2]<-sim.res[(ndat+1):(2*ndat),2]*(1-bias2)
var.calc<-NA
for(i in 1:nrow(sim.dat)){
	var.calc[i]<-var(sim.res[((i-1)*ndat+1):(i*ndat),2])
}
a <- lm(sim.res[,2]~sim.res[,1])
weights<-1/rep(var.calc,each=ndat)

rsq.bias5[g]<-summary(lm(sim.res[,2]~sim.res[,1],weights=weights))$r.squared
}
lines(density(rsq.bias5),col="red")





# make final combined plot for paper
####################################

bias.df <- data.frame(rsq=c(rsq,rsq.bias,rsq.bias2,rsq.bias4,rsq.bias5),bias=c(rep("No bias",10000),rep("Single bias 10%",10000),rep("Single bias 20%",10000),rep("Double bias 20% + 10%",10000),rep("Double bias 20% + 20%",10000)))
bias.df$bias<-factor(bias.df$bias, levels = levels(bias.df$bias)[c(3,4,5,1,2)])

library(ggplot2)
library(RColorBrewer)
ggplot(bias.df, aes(rsq, fill = bias, colour = bias)) +
  geom_density(alpha = 0.2, aes(linetype=bias), bw="SJ",trim=TRUE) +
  xlim(0.95, 1) +
  theme_minimal() +
  xlab(expression(R^2)) +
  theme(legend.title=element_blank()) +
  scale_colour_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") 

ggplot(bias.df, aes(y=rsq, x=factor(bias))) +
  geom_boxplot(alpha = 0.2) +
  ylim(0.98, 1) + coord_flip() +
  theme_minimal() +
  ylab(expression(R^2)) + xlab("") +
  theme(legend.position="none")



plot.list<-list(ggplot(bias.df1, aes(y=rsq, x=factor(bias))) +
  geom_boxplot(alpha = 0.2) +
  ylim(0.98, 1) + coord_flip() +
  theme_minimal() +
  ylab(expression(R^2)) + xlab("") +
  theme(legend.position="none"),
  
  ggplot(bias.df, aes(y=rsq, x=factor(bias))) +
  geom_boxplot(alpha = 0.2) +
  ylim(0.98, 1) + coord_flip() +
  theme_minimal() +
  ylab(expression(R^2)) + xlab("") +
  theme(legend.position="none")
)

plot_grid(plotlist=plot.list,nrow=2,labels=c("A","B"))