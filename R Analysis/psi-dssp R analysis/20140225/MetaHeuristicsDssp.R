#load some libraries ----------------------------------------------------------------

library(ggplot2)
library(GGally)
require(mgcv)
library(scales)
library(pastecs)
library(grid)
library(Hmisc)
library(MASS)
library(gridExtra)
library(EnvStats)
library(data.table)
library(reshape2)  

#end some libraries ----------------------------------------------------------------


#create some useful function--------------------------------------------------------

clc <- function() cat(rep("\n",50))  #clear console


quickStat <- function(x) {val<-c(stat.desc(x)[c("mean","std.dev","min","max")],skewness(x))
                          names(val)[5]<-"skewness"
                          print (val)
}

stderr <- function(x) sqrt(var(x)/length(x))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

bootStrapCorr <- function(x){
  R = 100                                    # the number of replicates
  cor1 <- numeric(R)
  cor2 <- numeric(R)      
  corr.values <- data.frame(cor1,cor2)                        # storage for the results
  for (i in 1:R) {
    group1 <- skinny[sample(1:nrow(skinny), 20000, replace=TRUE),]
    
    corr.values$cor1[i] <- cor.test(group1[,x],group1[,"pBestGap"])$estimate
    
    group2 <- skinny[sample(1:nrow(skinny), 20000, replace=TRUE),]
    corr.values$cor2[i] <- cor.test(group2[,x],group2$dsspGap)$estimate
  }
  t.test(corr.values$cor1,corr.values$cor2)
}


#end creation-----------------------------------------------------------------------

#change work directory--------------------------------------------------------------
getwd()
setwd("F:/OU/Research Projects/Active/Dynamic Slope Scaling Search/R Analysis/psi-dssp R analysis")
#end set new work directory---------------------------------------------------------

#************************************************************************************
#PSO analysis start
#************************************************************************************
#read data from csv file------------------------------------------------------------
PSOdata <- read.csv("PSOcsv.csv")  #PSO
head(PSOdata)
is.data.frame(PSOdata)
PSOdata <- data.table(PSOdata)

#end load data----------------------------------------------------------------------

#start PSO analysis------------------------------------------------------------
tables()  #check the table data list
PSOdata[,nodeCnt]#pick column nodeCnt

#general analysis ----------------------------------------------
#FCNF statical information-----------------
summary(PSOdata)
PSOdata[,c("density") := PSOdata[,noV/(2*choose(nodeCnt,2))]]  #density in all data setand creat new column
PSOdata[,density]
PSOdata[,c("psupplynodes") := PSOdata[,supplynodes/nodeCnt]]  #percentage of supply nodes
PSOdata[,c("pdemandnodes") := PSOdata[,demandnodes/nodeCnt]]  #percentage of demand nodes
PSOdata[,c("averagesupply") := PSOdata[,totalSupply/supplynodes]]   #caculate the average supply for each supply node
PSOdata[,c("gamma") := PSOdata[,(totalSupply/supplynodes)/avgFC]]   #gamma is defined in Dr.Nicholson's paper
mydata <- PSOdata
mydata$PSODSSPGap <- as.numeric(sub("%", "",mydata$PSODSSPGap))  #transfer str to number
is.numeric(mydata$PSODSSPGap)

#first try==============================================================================================================
#histogram of FCNF characteristics--
q1 <- ggplot(data=mydata,aes(x=density)) + geom_histogram(binwidth=0.1, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Density") +
  scale_y_continuous(labels = comma,limits=c(0,40)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

q2 <- ggplot(data=mydata,aes(x=psupplynodes)) + geom_histogram(binwidth=0.05, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = NULL,
       x = "Percentage of Supply Nodes") +
  scale_y_continuous(labels = comma,limits=c(0,80)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

q3 <- ggplot(data=mydata,aes(x=pdemandnodes)) + geom_histogram(binwidth=0.05, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Percentage of Demand Nodes") +
  scale_y_continuous(labels = comma,limits=c(0,80)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

q4 <- ggplot(data=mydata,aes(x=averagesupply)) + geom_histogram(binwidth=200, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = NULL,
       x = "Average Supply") +
  scale_y_continuous(labels = comma,limits=c(0,60)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

q5 <- ggplot(data=mydata,aes(x=avgFC)) + geom_histogram(binwidth=100, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Fixed to Variable Cost Ratio") +
  scale_y_continuous(labels = comma,limits=c(0,80)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

q6 <- ggplot(data=mydata,aes(x=gamma)) + geom_histogram(binwidth=.05, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = NULL,
       x = expression(paste(gamma," Characteristic"))) +
  scale_y_continuous(labels = comma,limits=c(0,60)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 


#save off histograms for paper -------------------------#
pdf(file="histogramSPhiRhoGamma.pdf",width=8,height=8)  #
grid.arrange(q1, q2, q3, q4, q5,q6,ncol=2)                    #  <------------------------------------------------------------------< save pdf images
dev.off()                                               #
#save off histograms for paper -------------------------#

#mean,std.dev,min,max,skewness--------------------------
quickStat(mydata$density)
quickStat(mydata$psupplynodes)
quickStat(mydata$pdemandnodes)
quickStat(mydata$averagesupply)
quickStat(mydata$avgFC)
quickStat(mydata$gamma)
mydata$PSODSSPGap <- as.numeric(sub("%", "",mydata$PSODSSPGap))  #transfer str to number
quickStat(mydata$PSODSSPGap)


#correlation test---------------------------------------
cor.test(mydata$density,mydata$PSODSSPGap,method="pearson")
cor.test(mydata$psupplynodes,mydata$PSODSSPGap,method="pearson")
cor.test(mydata$demandnodes,mydata$PSODSSPGap,method="pearson")
cor.test(mydata$averagesupply,mydata$PSODSSPGap,method="pearson")
cor.test(mydata$avgFC,mydata$PSODSSPGap,method="pearson")
cor.test(mydata$gamma,mydata$PSODSSPGap,method="pearson")
cor.test(mydata$PSODSSPGap,mydata$PSODSSPGap,method="pearson")

cor.test(mydata$nodeCnt,mydata$PSODSSPGap,method="pearson")
cor.test(mydata$noV,mydata$PSODSSPGap,method="pearson")

cor.test(mydata$density,mydata$PSObeatDSSP,method="pearson")
cor.test(mydata$psupplynodes,mydata$PSObeatDSSP,method="pearson")
cor.test(mydata$demandnodes,mydata$PSObeatDSSP,method="pearson")
cor.test(mydata$averagesupply,mydata$PSObeatDSSP,method="pearson")
cor.test(mydata$avgFC,mydata$PSObeatDSSP,method="pearson")
cor.test(mydata$gamma,mydata$PSObeatDSSP,method="pearson")
cor.test(mydata$PSODSSPGap,mydata$PSObeatDSSP,method="pearson")

cor.test(mydata$nodeCnt,mydata$PSObeatDSSP,method="pearson")
cor.test(mydata$noV,mydata$PSObeatDSSP,method="pearson")

#histogram for paper ----------------------------------------------------------------------- dsspGap histrogram
q1<-qplot(data=subset(mydata,PSODSSPGap>0),x=PSODSSPGap, binwidth=2) + theme_bw() +    
  labs(y = "Observations", x = "PSO-DSSP Gap (%)") +     
  geom_histogram(binwidth=2, color="black", fill = "light blue", size=.1) +        
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank()
  )  

#for paper! ----------------------------------#
pdf(file="dsspGap1.pdf",width=6,height=4)     #
print (q1)                                    # <------------------------------------------------------------------------< save pdf images
dev.off()                                     #
#for paper! ----------------------------------#

#PSO Performance Metrics and -----------------------------------------------------
#density------------------------------------------------------------------
#x<-quantile(mydata[,density],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("densityBins")] <- cut(mydata[,density], unique(quantile(mydata[,density],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,density, densityBins, PSODSSPGap )])
is.numeric(dt$PSODSSPGap)

dt2<-data.frame(dt[,list(obs=length(seed),meandensity=mean(density),meanGap=mean(PSODSSPGap),sdGap=sd(PSODSSPGap),seGap=stderr(PSODSSPGap)),by=densityBins])
#dt2<-data.frame()
dt2$meandensity

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g1<-ggplot(data=dt2,aes(x=meandensity,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,15)) +
  labs(y = "Mean PSO vs DSSP Gap (%)", x = "Average Density") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")


#percentage of supply nodes----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("psupplynodesBins")] <- cut(mydata[,psupplynodes], unique(quantile(mydata[,psupplynodes],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,psupplynodes,psupplynodesBins, PSODSSPGap )])
is.numeric(dt$PSODSSPGap)

dt2<-data.frame(dt[,list(obs=length(seed),meanpsupplynodes=mean(psupplynodes),meanGap=mean(PSODSSPGap),sdGap=sd(PSODSSPGap),seGap=stderr(PSODSSPGap)),by=psupplynodesBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g2<-ggplot(data=dt2,aes(x=meanpsupplynodes,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,15)) +
  labs(y ="Mean PSO vs DSSP Gap (%)" , x = "Percentage of Supply Nodes") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")


#percentage of demand nodes----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("pdemandnodesBins")] <- cut(mydata[,pdemandnodes], unique(quantile(mydata[,pdemandnodes],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,pdemandnodes,pdemandnodesBins, PSODSSPGap )])
is.numeric(dt$PSODSSPGap)

dt2<-data.frame(dt[,list(obs=length(seed),meanpdemandnodes=mean(pdemandnodes),meanGap=mean(PSODSSPGap),sdGap=sd(PSODSSPGap),seGap=stderr(PSODSSPGap)),by=pdemandnodesBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g3<-ggplot(data=dt2,aes(x=meanpdemandnodes,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,15)) +
  labs(y ="Mean PSO vs DSSP Gap (%)" , x = "Percentage of Demand Nodes") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")


#average supply----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("averagesupplyBins")] <- cut(mydata[,averagesupply], unique(quantile(mydata[,averagesupply],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,averagesupply,averagesupplyBins, PSODSSPGap )])
is.numeric(dt$PSODSSPGap)

dt2<-data.frame(dt[,list(obs=length(seed),meanaveragesupply=mean(averagesupply),meanGap=mean(PSODSSPGap),sdGap=sd(PSODSSPGap),seGap=stderr(PSODSSPGap)),by=averagesupplyBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g4<-ggplot(data=dt2,aes(x=meanaveragesupply,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,15)) +
  labs(y ="Mean PSO vs DSSP Gap (%)", x = "Average Supply") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

#F/C ratio----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("avgFCBins")] <- cut(mydata[,avgFC], unique(quantile(mydata[,avgFC],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,avgFC,avgFCBins, PSODSSPGap )])
is.numeric(dt$PSODSSPGap)

dt2<-data.frame(dt[,list(obs=length(seed),meanavgFC = mean(avgFC),meanGap=mean(PSODSSPGap),sdGap=sd(PSODSSPGap),seGap=stderr(PSODSSPGap)),by=avgFCBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g5<-ggplot(data=dt2,aes(x=meanavgFC,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,18)) +
  labs(y ="Mean PSO vs DSSP Gap (%)", x = "Fixed to Variable Cost Ratio") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

#gamma----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("gammaBins")] <- cut(mydata[,gamma], unique(quantile(mydata[,gamma],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,gamma,gammaBins, PSODSSPGap )])
is.numeric(dt$PSODSSPGap)

dt2<-data.frame(dt[,list(obs=length(seed),meangamma = mean(gamma),meanGap=mean(PSODSSPGap),sdGap=sd(PSODSSPGap),seGap=stderr(PSODSSPGap)),by=gammaBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g6<-ggplot(data=dt2,aes(x=meangamma,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,18)) +
  labs(y ="Mean PSO vs DSSP Gap (%)", x =expression(paste(gamma," Characteristic"))) +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

#for paper! ----------------------------------------#
pdf(file="PSOGapPerformance.pdf",width=8,height=8)     #
grid.arrange(g1,g2,g3,g4,g5,g6,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

#PSO Performance Metrics (percentage of PSOwin) and -----------------------------------------------------
#density------------------------------------------------------------------
#x<-quantile(mydata[,density],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("densityBins")] <- cut(mydata[,density], unique(quantile(mydata[,density],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,density, densityBins, PSObeatDSSP )])

dt2<-data.frame(dt[,list(obs=length(seed),meandensity=mean(density),meanGap=mean(PSObeatDSSP),sdGap=sd(PSObeatDSSP),seGap=stderr(PSObeatDSSP)),by=densityBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

gg1<-ggplot(data=dt2,aes(x=meandensity,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0.6,1)) +
  labs(y = "Percentage of PSO beats DSSP", x = "Average Density") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

#percentage of supply nodes----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("psupplynodesBins")] <- cut(mydata[,psupplynodes], unique(quantile(mydata[,psupplynodes],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,psupplynodes,psupplynodesBins, PSObeatDSSP)])

dt2<-data.frame(dt[,list(obs=length(seed),meanpsupplynodes=mean(psupplynodes),meanGap=mean(PSObeatDSSP),sdGap=sd(PSObeatDSSP),seGap=stderr(PSObeatDSSP)),by=psupplynodesBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

gg2<-ggplot(data=dt2,aes(x=meanpsupplynodes,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0.6,1.0)) +
  labs(y ="Percentage of PSO beats DSSP" , x = "Percentage of Supply Nodes") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")


#percentage of demand nodes----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("pdemandnodesBins")] <- cut(mydata[,pdemandnodes], unique(quantile(mydata[,pdemandnodes],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,pdemandnodes,pdemandnodesBins, PSObeatDSSP )])

dt2<-data.frame(dt[,list(obs=length(seed),meanpdemandnodes=mean(pdemandnodes),meanGap=mean(PSObeatDSSP),sdGap=sd(PSObeatDSSP),seGap=stderr(PSObeatDSSP)),by=pdemandnodesBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

gg3<-ggplot(data=dt2,aes(x=meanpdemandnodes,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0.6,1.0)) +
  labs(y ="Percentage of PSO beats DSSP" , x = "Percentage of Demand Nodes") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")


#average supply----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("averagesupplyBins")] <- cut(mydata[,averagesupply], unique(quantile(mydata[,averagesupply],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,averagesupply,averagesupplyBins, PSObeatDSSP )])
is.numeric(dt$PSObeatDSSP)

dt2<-data.frame(dt[,list(obs=length(seed),meanaveragesupply=mean(averagesupply),meanGap=mean(PSObeatDSSP),sdGap=sd(PSObeatDSSP),seGap=stderr(PSObeatDSSP)),by=averagesupplyBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

gg4<-ggplot(data=dt2,aes(x=meanaveragesupply,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0.6,1.0)) +
  labs(y ="Percentage of PSO beats DSSP", x = "Average Supply") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

#F/C ratio----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("avgFCBins")] <- cut(mydata[,avgFC], unique(quantile(mydata[,avgFC],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,avgFC,avgFCBins, PSObeatDSSP )])


dt2<-data.frame(dt[,list(obs=length(seed),meanavgFC = mean(avgFC),meanGap=mean(PSObeatDSSP),sdGap=sd(PSObeatDSSP),seGap=stderr(PSObeatDSSP)),by=avgFCBins])
dt2
limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

gg5<-ggplot(data=dt2,aes(x=meanavgFC,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0.6,1.0)) +
  labs(y ="Percentage of PSO beats DSSP", x = "Fixed to Variable Cost Ratio") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

#gamma----------------------------------------------------------
#x<-quantile(mydata[,psupplynodes],  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
#x["0%"]<- x["0%"]*.9
mydata[,c("gammaBins")] <- cut(mydata[,gamma], unique(quantile(mydata[,gamma],probs=seq(0,1,0.1))),include.lowest=TRUE)

dt<-data.table(mydata[,list(seed,gamma,gammaBins, PSObeatDSSP )])
is.numeric(dt$PSObeatDSSP)

dt2<-data.frame(dt[,list(obs=length(seed),meangamma = mean(gamma),meanGap=mean(PSObeatDSSP),sdGap=sd(PSObeatDSSP),seGap=stderr(PSObeatDSSP)),by=gammaBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

gg6<-ggplot(data=dt2,aes(x=meangamma,y=meanGap)) + 
  geom_point(size=2,shape=6) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0.6,1.0)) +
  labs(y ="Percentage of PSO beats DSSP", x =expression(paste(gamma," Characteristic"))) +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

#for paper! ----------------------------------------#
pdf(file="PSOWinPerformance.pdf",width=8,height=8)     #
grid.arrange(gg1,gg2,gg3,gg4,gg5,gg6,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

#for paper! ----------------------------------------#
pdf(file="PSOWinGapPerformance1.pdf",width=8,height=8)     #
grid.arrange(g1,gg1,g2,gg2,g3,gg3,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

#for paper! ----------------------------------------#
pdf(file="PSOWinGapPerformance2.pdf",width=8,height=8)     #
grid.arrange(g4,gg4,g5,gg5,g6,gg6,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

#boxplot GAP and win based on nodes and noV-----------------------------------------------------------
#for paper! ----------------------------------------#
pdf(file="PSOBoxplot.pdf",width=6,height=4)     #
boxplot(PSODSSPGap~nodeCnt,data=mydata, main="PSO vs DSSP Gap Boxplot", 
        xlab="Number of Nodes", ylab="Gap(%)")
# <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

#PSO Performance and noV-----------------------------------------------------
pdf(file="scatterplotnov .pdf",width=6,height=4) 
qplot(noV,PSODSSPGap,data=mydata, facets=NULL, main="PSO vs DSSP Gap", xlab="Number of Variables",ylab="Gap (%)")
dev.off()

#psi value study--------------------------------------------------------------------
#histogram for paper ----------------------------------------------------------------------- dsspGap histrogram
q1<-qplot(data=subset(mydata,PSODSSPGap>0),x=PSOpsi, binwidth=0.1) + theme_bw() +    
  labs(y = "Observations", x = expression(paste(psi, " value"))) + 
  scale_y_continuous(labels = comma, limits=c(0,30)) +
  geom_histogram(binwidth=0.1, color="black", fill = "light blue", size=.1) +        
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank()
  )  

#for paper! ----------------------------------#
pdf(file="psihis.pdf",width=6,height=4)     #
print (q1)                                    # <------------------------------------------------------------------------< save pdf images
dev.off()                                     #
#for paper! ----------------------------------#

#correlation analysis
cor.test(mydata$density,mydata$PSOpsi,method="pearson")
cor.test(mydata$psupplynodes,mydata$PSOpsi,method="pearson")
cor.test(mydata$demandnodes,mydata$PSOpsi,method="pearson")
cor.test(mydata$averagesupply,mydata$PSOpsi,method="pearson")
cor.test(mydata$avgFC,mydata$PSOpsi,method="pearson")
cor.test(mydata$gamma,mydata$PSOpsi,method="pearson")
cor.test(mydata$PSODSSPGap,mydata$PSOpsi,method="pearson")
cor.test(mydata$nodeCnt,mydata$PSOpsi,method="pearson")
cor.test(mydata$noV,mydata$PSOpsi,method="pearson")
#end psi study----------------------------------------------------------------------------

#end PSO analysis--------------------------------------------------------
#end first try==============================================================================================================

#second try===========================================================================================
#scatter plot with smooth curve
mydata$PSODSSPGap <- as.numeric(sub("%", "",mydata$PSODSSPGap))  #transfer str to number
quickStat(mydata$PSODSSPGap)
is.numeric(mydata$PSODSSPGap)
p1 <- qplot(density, PSODSSPGap, data = mydata,geom = c("point", "smooth"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Density")
p11 <- p1 + facet_wrap(~nodeCnt)

p2 <- qplot(psupplynodes, PSODSSPGap, data = mydata,geom = c("point", "smooth"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Percentage of Supply nodes")
p22 <- p2 + facet_wrap(~nodeCnt)
p3 <- qplot(pdemandnodes, PSODSSPGap, data = mydata,geom = c("point", "smooth"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Percentage of Supply nodes")
p33 <- p3 + facet_wrap(~nodeCnt)
#for paper! ----------------------------------------#
pdf(file="ScatterPlotSmoothCurve1.pdf",width=8,height=8)     #
grid.arrange(p1,p11,p2,p22,p3,p33,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

p1 <- qplot(averagesupply, PSODSSPGap, data = mydata,geom = c("point", "smooth"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Average Supply")
p11 <- p1 + facet_wrap(~nodeCnt)

p2 <- qplot(avgFC, PSODSSPGap, data = mydata,geom = c("point", "smooth"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Fixed to Variable Cost Ratio")
p22 <- p2 + facet_wrap(~nodeCnt)
p3 <- qplot(gamma, PSODSSPGap, data = mydata,geom = c("point", "smooth"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x =  expression(paste(gamma," Characteristic")))
p33 <- p3 + facet_wrap(~nodeCnt)
#for paper! ----------------------------------------#
pdf(file="ScatterPlotSmoothCurve2.pdf",width=8,height=8)     #
grid.arrange(p1,p11,p2,p22,p3,p33,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

p1 <- qplot(noV, PSODSSPGap, data = mydata,geom = c("point", "smooth"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x =  "Number of Variables")
p11 <- p1 + facet_wrap(~nodeCnt)

#for paper! ----------------------------------------#
pdf(file="ScatterPlotSmoothCurve3.pdf",width=8,height=6)     #
grid.arrange(p11)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

#correlation analysis
cor.test(mydata$density,mydata$PSODSSPGap,method="spearm")
cor.test(mydata$psupplynodes,mydata$PSODSSPGap,method="spearm")
cor.test(mydata$demandnodes,mydata$PSODSSPGap,method="spearm")
cor.test(mydata$averagesupply,mydata$PSODSSPGap,method="spearm")
cor.test(mydata$avgFC,mydata$PSODSSPGap,method="spearm")
cor.test(mydata$gamma,mydata$PSODSSPGap,method="spearm")

#histogram of PSO beat DSSP based on number of nodes
p1 <- ggplot(data=mydata,aes(x=PSObeatDSSP)) + geom_histogram(binwidth=0.1, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "PSO overperform DSSP") +
  scale_y_continuous(labels = comma,limits=c(0,70)) +
  scale_x_continuous(labels = comma,limits=c(0,1)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

p11 <- p1+ facet_wrap(~nodeCnt)

#for paper! ----------------------------------------#
pdf(file="HistofPSObeatDSSP.pdf",width=8,height=6)     #
grid.arrange(p11)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

#logistic regression model analyze the possibility of PSO overperfom DSSP---------------
sink("logistic regression with nodeCnt.txt")
glm.out <- glm(PSObeatDSSP ~ nodeCnt, family=binomial(logit), data=mydata)
summary(glm.out)
sink()

sink("logistic regression with all variables.txt")
glm.out <- glm(PSObeatDSSP ~ nodeCnt+density+noV+avgFC+psupplynodes+pdemandnodes+averagesupply+gamma, family=binomial(logit), data=mydata)
summary(glm.out)
glm.out <- glm(PSObeatDSSP ~ nodeCnt*density*noV*avgFC*psupplynodes*pdemandnodes*averagesupply*gamma, family=binomial(logit), data=mydata)
summary(glm.out)
sink()

sink("logistic regression nodeCnt=25.txt")
mydata1 <- mydata[nodeCnt==25,]
glm.out <- glm(PSObeatDSSP ~ density+noV+avgFC+psupplynodes+pdemandnodes+averagesupply+gamma, family=binomial(logit), data=mydata)
summary(glm.out)
glm.out <- glm(PSObeatDSSP ~ density*noV*avgFC*psupplynodes*pdemandnodes*averagesupply*gamma, family=binomial(logit), data=mydata1)
summary(glm.out)
sink()

sink("logistic regression nodeCnt=50.txt")
mydata1 <- mydata[nodeCnt==50,]
glm.out <- glm(PSObeatDSSP ~ density+noV+avgFC+psupplynodes+pdemandnodes+averagesupply+gamma, family=binomial(logit), data=mydata)
summary(glm.out)
glm.out <- glm(PSObeatDSSP ~ density*noV*avgFC*psupplynodes*pdemandnodes*averagesupply*gamma, family=binomial(logit), data=mydata1)
summary(glm.out)
sink()

sink("logistic regression nodeCnt=100.txt")
mydata1 <- mydata[nodeCnt==100,]
glm.out <- glm(PSObeatDSSP ~ density+noV+avgFC+psupplynodes+pdemandnodes+averagesupply+gamma, family=binomial(logit), data=mydata)
summary(glm.out)
glm.out <- glm(PSObeatDSSP ~ density*noV*avgFC*psupplynodes*pdemandnodes*averagesupply*gamma, family=binomial(logit), data=mydata1)
summary(glm.out)
sink()
#end second try===========================================================================================


#Third try===========================================================================================
#step-wise ------------------------------------------------------------------------------
sink("step-wise.txt")
null <- glm(PSObeatDSSP ~ 1, family=binomial(logit), data=mydata)
full <- glm(PSObeatDSSP ~ nodeCnt*density*noV*avgFC*psupplynodes*pdemandnodes*averagesupply*gamma, family=binomial(logit), data=mydata)
glm.step <- stepAIC(null, scope=list(lower=null, upper=full), direction="both")
summary(glm.step)
glm.step$anova
sink()

#time------------------------------------------------------------------------------------
head(mydata$PSOtime)
b1 <- boxplot(PSOtime~nodeCnt, data=mydata[nodeCnt==25,], main="PSO Time Boxplot", 
              xlab="Number of Nodes", ylab="Time(s)")
b2 <- boxplot(PSOtime~nodeCnt, data=mydata[nodeCnt==50,], main="PSO Time Boxplot", 
              xlab="Number of Nodes", ylab="Time(s)")
b3 <- boxplot(PSOtime~nodeCnt, data=mydata[nodeCnt==100,], main="PSO Time Boxplot", 
              xlab="Number of Nodes", ylab="Time(s)")
# <------------------------------------------------------------------------< save pdf images
pdf(file="PSOTimeBoxplot.pdf",width=8,height=4)     #
par(mfrow=c(1,3))
boxplot(PSOtime~nodeCnt, data=mydata[nodeCnt==25,], main="", 
        xlab="Nodes =25", ylab="Time(s)")
boxplot(PSOtime~nodeCnt, data=mydata[nodeCnt==50,], main="", 
        xlab="Nodes=50", ylab="")
boxplot(PSOtime~nodeCnt, data=mydata[nodeCnt==100,], main="", 
        xlab="Nodes=100", ylab="")
dev.off()                                           #
par(mfrow=c(1,1))
#for paper! ----------------------------------------#
#narow the scope where PSO has a relative high performance-------------------------------
mydata$PSODSSPGap <- as.numeric(sub("%", "",mydata$PSODSSPGap))  #transfer str to number
is.numeric(mydata$PSODSSPGap)
summary(mydata$PSODSSPGap)
mean(mydata$PSODSSPGap)

mydata2 <- subset(mydata, PSODSSPGap>15)
head(mydata2$density)
p1 <- qplot(density, PSODSSPGap, data = mydata2,geom = c("point"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Density")
p11 <- p1 + facet_wrap(~nodeCnt)

p2 <- qplot(psupplynodes, PSODSSPGap, data = mydata2,geom = c("point"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Percentage of Supply nodes")
p22 <- p2 + facet_wrap(~nodeCnt)
p3 <- qplot(pdemandnodes, PSODSSPGap, data = mydata2,geom = c("point"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Percentage of Demand nodes")
p33 <- p3 + facet_wrap(~nodeCnt)
#for paper! ----------------------------------------#
pdf(file="ScatterPlot1highgap.pdf",width=8,height=8)     #
grid.arrange(p1,p11,p2,p22,p3,p33,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

p1 <- qplot(averagesupply, PSODSSPGap, data = mydata2,geom = c("point"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Average Supply")
p11 <- p1 + facet_wrap(~nodeCnt)

p2 <- qplot(avgFC, PSODSSPGap, data = mydata2,geom = c("point"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x = "Fixed to Variable Cost Ratio")
p22 <- p2 + facet_wrap(~nodeCnt)
p3 <- qplot(gamma, PSODSSPGap, data = mydata2,geom = c("point"))+ theme_bw() +    
  labs(y = "PSO vs DSSP Gap(%)", x =  expression(paste(gamma," Characteristic")))
p33 <- p3 + facet_wrap(~nodeCnt)
#for paper! ----------------------------------------#
pdf(file="ScatterPlot2highgap.pdf",width=8,height=8)     #
grid.arrange(p1,p11,p2,p22,p3,p33,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

#histgram---------------------------------
p1 <- ggplot(data=mydata2,aes(x=density)) + geom_histogram(binwidth=0.2, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Density") +
  #scale_y_continuous(labels = comma,limits=c(0,18)) +
  scale_x_continuous(labels = comma,limits=c(0,1)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

p11 <- p1+ facet_wrap(~nodeCnt)

p2 <- ggplot(data=mydata2,aes(x=psupplynodes)) + geom_histogram(binwidth=0.02, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Percentage of Supply Nodes") +
  #scale_y_continuous(labels = comma,limits=c(0,18)) +
  scale_x_continuous(labels = comma,limits=c(0,0.7)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

p22 <- p2+ facet_wrap(~nodeCnt)

p3 <- ggplot(data=mydata2,aes(x=pdemandnodes)) + geom_histogram(binwidth=0.02, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Percentage of Demand Nodes") +
  #scale_y_continuous(labels = comma,limits=c(0,18)) +
  scale_x_continuous(labels = comma,limits=c(0,0.7)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

p33 <- p3 + facet_wrap(~nodeCnt)

#for paper! ----------------------------------------#
pdf(file="hist1highgap.pdf",width=8,height=8)     #
grid.arrange(p1,p11,p2,p22,p3,p33,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

p1 <- ggplot(data=mydata2,aes(x=averagesupply)) + geom_histogram(binwidth=100, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Average Supply") +
  #scale_y_continuous(labels = comma,limits=c(0,18)) +
  #scale_x_continuous(labels = comma,limits=c(0,1)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

p11 <- p1+ facet_wrap(~nodeCnt)

p2 <- ggplot(data=mydata2,aes(x=avgFC)) + geom_histogram(binwidth=50, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "F/C Ratio") +
  #scale_y_continuous(labels = comma,limits=c(0,18)) +
  #scale_x_continuous(labels = comma,limits=c(0,1)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

p22 <- p2+ facet_wrap(~nodeCnt)

p3 <- ggplot(data=mydata2,aes(x=gamma)) + geom_histogram(binwidth=0.05, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Gamma") +
  #scale_y_continuous(labels = comma,limits=c(0,18)) +
  scale_x_continuous(labels = comma,limits=c(0,0.7)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

p33 <- p3 + facet_wrap(~nodeCnt)

#for paper! ----------------------------------------#
pdf(file="hist2highgap.pdf",width=8,height=8)     #
grid.arrange(p1,p11,p2,p22,p3,p33,ncol=2)              # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#

#double histogram
p1 <- ggplot(data=mydata2,aes(x=density)) + geom_histogram(binwidth=0.1, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Density") +
  #scale_y_continuous(labels = comma,limits=c(0,18)) +
  scale_x_continuous(labels = comma,limits=c(0,1)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

p11 <- p1+ facet_wrap(~nodeCnt)

p2 <- ggplot(data=mydata,aes(x=density)) + geom_histogram(binwidth=0.1, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
  labs(y = "Observations",
       x = "Density") +
  #scale_y_continuous(labels = comma,limits=c(0,18)) +
  scale_x_continuous(labels = comma,limits=c(0,1)) +
  theme(axis.ticks = element_blank(), 
        # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

p11 <- p1+ facet_wrap(~nodeCnt)

par(mfrow=c(1,1))

p1<-hist(mydata$density)
p2<-hist(mydata2$density)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1), main="Histograme of Density",xlab="Density")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)  # second

p1<-hist(mydata$avgFC)
p2<-hist(mydata2$avgFC)
plot( p1, col=rgb(0,0,1,1/4), main="Histograme of F/C Ratio",xlab="Fixed vs Variable Ratio")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), add=T)  # second

p1<-hist(mydata$psupplynodes)
p2<-hist(mydata2$psupplynodes)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1), main="Histograme of % Supply Nodes",xlab="Supply Nodes")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)  # second

p1<-hist(mydata$pdemandnodes)
p2<-hist(mydata2$pdemandnodes)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1), main="Histograme of % Demand Nodes",xlab="Demand Nodes")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)  # second

p1<-hist(mydata$averagesupply)
p2<-hist(mydata2$averagesupply)
plot( p1, col=rgb(0,0,1,1/4), main="Histograme of Average Supply",xlab="Average Supply")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), add=T)  # second

p1<-hist(mydata$gamma)
p2<-hist(mydata2$gamma)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1), main="Histograme of Gamma",xlab="Gamma")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)  # second

#end third try===========================================================================================

#************************************************************************************
#PSO analysis end
#************************************************************************************


#************************************************************************************
#SA-B analysis start
#************************************************************************************
#read data from csv file------------------------------------------------------------
SABdata <- read.csv("SABcsv.csv")  #SAB
head(SABdata)
is.data.frame(SABdata)
SABdata <- data.table(SABdata)

#end load data----------------------------------------------------------------------

tables()  #check the table data list
SABdata[,nodeCnt]#pick column nodeCnt

#general analysis ----------------------------------------------
#FCNF statical information-----------------
summary(SABdata)
SABdata[,c("density") := SABdata[,noV/(2*choose(nodeCnt,2))]]  #density in all data setand creat new column
SABdata[,density]
SABdata[,c("psupplynodes") := SABdata[,supplynodes/nodeCnt]]  #percentage of supply nodes
SABdata[,c("pdemandnodes") := SABdata[,demandnodes/nodeCnt]]  #percentage of demand nodes
SABdata[,c("averagesupply") := SABdata[,totalSupply/supplynodes]]   #caculate the average supply for each supply node
SABdata[,c("gamma") := SABdata[,(totalSupply/supplynodes)/avgFC]]   #gamma is defined in Dr.Nicholson's paper
mydata <- SABdata
mydata$SABDSSPGap <- as.numeric(sub("%", "",mydata$SABDSSPGap))  #transfer str to number
is.numeric(mydata$SABDSSPGap)
head(mydata)
boxplot(SABDSSPGap~nodeCnt,data=mydata, main="SA-B vs DSSP Gap Boxplot", 
        xlab="Number of Nodes", ylab="Gap(%)")
par(mfrow=c(1,3))
boxplot(SABtime~nodeCnt, data=mydata[nodeCnt==25,], main="", 
        xlab="Nodes =25", ylab="Time(s)")
boxplot(SABtime~nodeCnt, data=mydata[nodeCnt==50,], main="", 
        xlab="Nodes=50", ylab="")
boxplot(SABtime~nodeCnt, data=mydata[nodeCnt==100,], main="", 
        xlab="Nodes=100", ylab="")
par(mfrow=c(1,1))




#************************************************************************************
#SA-VF analysis start
#************************************************************************************
SAVFdata <- read.csv("SAVFcsv.csv")  #SAB
head(SAVFdata)
is.data.frame(SAVFdata)
SAVFdata <- data.table(SAVFdata)

#end load data----------------------------------------------------------------------

tables()  #check the table data list
SAVFdata[,nodeCnt]#pick column nodeCnt

#general analysis ----------------------------------------------
#FCNF statical information-----------------
summary(SAVFdata)
SAVFdata[,c("density") := SAVFdata[,noV/(2*choose(nodeCnt,2))]]  #density in all data setand creat new column
SAVFdata[,density]
SAVFdata[,c("psupplynodes") := SAVFdata[,supplynodes/nodeCnt]]  #percentage of supply nodes
SAVFdata[,c("pdemandnodes") := SAVFdata[,demandnodes/nodeCnt]]  #percentage of demand nodes
SAVFdata[,c("averagesupply") := SAVFdata[,totalSupply/supplynodes]]   #caculate the average supply for each supply node
SAVFdata[,c("gamma") := SAVFdata[,(totalSupply/supplynodes)/avgFC]]   #gamma is defined in Dr.Nicholson's paper
mydata <- SAVFdata
mydata$SAVFDSSPGap <- as.numeric(sub("%", "",mydata$SAVFDSSPGap))  #transfer str to number
is.numeric(mydata$SAVFDSSPGap)
head(mydata)
boxplot(SAVFDSSPGap~nodeCnt,data=mydata, main="SAVF vs DSSP Gap Boxplot", 
        xlab="Number of Nodes", ylab="Gap(%)")
par(mfrow=c(1,3))
boxplot(SAVFtime~nodeCnt, data=mydata[nodeCnt==25,], main="", 
        xlab="Nodes =25", ylab="Time(s)")
boxplot(SAVFtime~nodeCnt, data=mydata[nodeCnt==50,], main="", 
        xlab="Nodes=50", ylab="")
boxplot(SAVFtime~nodeCnt, data=mydata[nodeCnt==100,], main="", 
        xlab="Nodes=100", ylab="")
par(mfrow=c(1,1))

#ANOVA for gap and time based on PSO SAB SAVF
#load new data set
mydata <- read.csv("GapTimecsv.csv")
head(mydata)
mydata <- data.table(mydata)
is.data.table(mydata)
Algorithm <- factor(mydata$Algorithm)
#separate with nodes
data1 <- mydata[mydata$Node==25,]
data2 <- mydata[mydata$Node==50,]
data3 <- mydata[mydata$Node==100,]

sink("ANOVA and Tukey.txt")
print("Node=25")
print("ANOVA of Gap")
aov.out <- aov(Gap ~ Algorithm , data=data1)
summary(aov.out)
print("ANOVA of Time")
aov.out <- aov(Time ~ Algorithm , data=data1)
summary(aov.out)
print("TukeyTest")
TukeyHSD(aov.out)

print("Node=50")
print("ANOVA of Gap")
aov.out <- aov(Gap ~ Algorithm , data=data2)
summary(aov.out)
print("ANOVA of Time")
aov.out <- aov(Time ~ Algorithm , data=data2)
summary(aov.out)
print("TukeyTest")
TukeyHSD(aov.out)

print("Node=100")
print("ANOVA of Gap")
aov.out <- aov(Gap ~ Algorithm , data=data3)
summary(aov.out)
print("ANOVA of Time")
aov.out <- aov(Time ~ Algorithm , data=data3)
summary(aov.out)
print("TukeyTest")
TukeyHSD(aov.out)
sink()