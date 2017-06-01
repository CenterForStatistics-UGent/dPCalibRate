library(reshape2)

load('power.RData')
plot.suppl<-list()
par(mfrow=c(1,5))
typeI<-data.frame(OLS=power.noweights[,1],WLS=power[,1],WLS.HCCM=power.hccm[,1])
typeI<-melt(typeI)
typeI<-data.frame(typeI,replicates=c(rep(2:8,3)))
colnames(typeI)<-c("type","size","replicates")
plot.suppl[[1]]<-ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05) + ggtitle("Quadratic effect test") +theme(legend.position="none")


load('powerruns.RData')

typeI<-data.frame(OLS=power.noweights[,1],WLS=power[,1],WLS.HCCM=power.hccm[,1])
typeI<-melt(typeI)
typeI<-data.frame(typeI,replicates=c(rep(2:8,3)))
colnames(typeI)<-c("type","size","replicates")
plot.suppl[[2]]<-ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05)+ ggtitle("Runs test") +theme(legend.position="none")


load('powerfreq.RData')

typeI<-data.frame(OLS=power.noweights[,1],WLS=power[,1],WLS.HCCM=power.hccm[,1])
typeI<-melt(typeI)
typeI<-data.frame(typeI,replicates=c(rep(2:8,3)))
colnames(typeI)<-c("type","size","replicates")
plot.suppl[[3]]<-ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05)+ ggtitle("Frequency within a\nblock test") +theme(legend.position="none")


load('lackoffit.RData')

power.hccm<-power
typeI<-data.frame(OLS=power.noweights[,1],WLS=power[,1],WLS.HCCM=power.hccm[,1])
typeI<-melt(typeI)
typeI<-data.frame(typeI,replicates=c(rep(2:8,3)))
colnames(typeI)<-c("type","size","replicates")
plot.suppl[[4]]<-ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05)+ ggtitle("Lack of fit test") +theme(legend.position="none")


load('poweraccuracy.RData')

typeI<-data.frame(OLS=power.noweights[,1],WLS=power[,1],WLS.HCCM=power.hccm[,1])
typeI<-melt(typeI)
typeI<-data.frame(typeI,replicates=c(rep(2:8,3)))
colnames(typeI)<-c("type","size","replicates")
plot.suppl[[5]]<-ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05)+ ggtitle("Accuracy test") +theme(legend.position="none")

#get legend only
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}
plot.temp<-ggplot(typeI,aes(x=replicates,y=size,group=type,colour=type)) + geom_line(aes(linetype=type)) + theme_minimal() + geom_hline(yintercept=0.05)+ ggtitle("Accuracy test")

plot.suppl[[6]] <- g_legend(plot.temp)

plot_grid(plotlist=plot.suppl, align='h', labels=c('A', 'B', 'C', 'D','E',''), nrow=2)