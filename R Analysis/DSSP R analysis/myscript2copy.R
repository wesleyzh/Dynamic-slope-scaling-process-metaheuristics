


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

#create some useful functions --------------------------------------------------------

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
stderr <- function(x) sqrt(var(x)/length(x))

quickStat <- function(x) {val<-c(stat.desc(x)[c("mean","std.dev","min","max")],skewness(x))
                          names(val)[5]<-"skewness"
                          print (val)
                          }

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

#read in data; merge and etc. as well ------------------------------------------------

mydata <- read.csv("~/Research and Grants/Research/Active/Time Varying Cost/DSSP search results/Phi_Density_Scan_12 -FINAL 131912.txt")

metadata <- read.csv("~/Research and Grants/Research/Active/Time Varying Cost/DSSP search results/metadata40000.txt")
metadata$noV<-NULL
metadata$rhsMin<-NULL
metadata$rhsMax<-NULL
metadata$cMin<-NULL
metadata$cMax<-NULL
metadata$fMin<-NULL
metadata$fMax<-NULL

mydata <- merge(mydata,metadata,by="seed")

scan2 <- read.csv("~/Research and Grants/Research/Time Varying Cost/DSSP search results/Phi_Density_Scan_12 from 2.txt")
scan2$pBest2<-scan2$pBest
scan2$pBestIt2<-scan2$pBestIt
scan2$pBestTime2<-scan2$pBestTime
scan2<-scan2[,c("seed","pBest2","pBestIt2","pBestTime2")]

mydata <- merge(mydata,scan2,by="seed")

#create primary variables for analysis -----------------------------------------------------

mydata$rhoS<-(mydata$totSupply/mydata$supplyNodes)
mydata$rhoD<-(mydata$totSupply/mydata$demandNodes)
mydata$rho<- (mydata$rhoS+ mydata$rhoD)/2

summary(mydata$rho)

#mydata$ravg<-(mydata$rhsMin+mydata$rhsMax)/2
#mydata$rhoMax<-mapply(FUN=max, mydata$rhoS, mydata$rhoD, na.rm=TRUE)

mydata$density = (mydata$noV/(12*11))*100

mydata$dsspGap = 100*(mydata$dsspObj-mydata$optObj)/mydata$optObj
mydata$dsspGapI <- (mydata$dsspGap>0.0001)

mydata$pBestGap = 100*(mydata$pBestObj-mydata$optObj)/mydata$optObj
mydata$pBestGapI <- (mydata$pBestGap>0)
mydata$pBestGapI <- (mydata$pBestGap>0.0001)

mydata$pBestIL <- (mydata$pBest < 1)
mydata$pBest2IL <- (mydata$pBest2 < 1)

mydata$pBestBetterI <- (mydata$pBestGap < mydata$dsspGap)

mydata$Improvement <- mydata$dsspGap - mydata$pBestGap

mydata$phi<-mydata$avgFC
mydata$avgFC <- NULL
#mydata$oldphi = (mydata$fMax - mydata$fMin)/(mydata$cMax - mydata$cMin)

mydata$gamma= mydata$rho/mydata$phi
mydata$logGamma = log10(mydata$gamma)


mydata$phiBin<-cut(mydata$phi, breaks= c(0,1000,max(mydata$phi)+1), labels = F)
                   )






#histograms for paper -----------------------------------------------------------------------
q1<-ggplot(data=mydata,aes(x=totSupply)) +
  geom_histogram(binwidth=450, color="black", fill = "light blue", size=.1) +
  #geom_density(alpha=0.1, fill="#FF6666") +
  # geom_vline(aes(xintercept=mean(totSupply, na.rm=T)),  color="black", linetype="dashed", size=1) +
  labs(y = "Observations",
       x = "Total Supply") + 
       #x = expression(paste("Total Supply (", italic(S),")"))) + 
  scale_y_continuous(labels = comma, limits=c(0,6000)) +
  scale_x_continuous(labels = comma, limits=c(0,12000)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2))  

q2<-ggplot(data=mydata,aes(x=phi)) + 
  geom_histogram(binwidth=500, color="black", fill = "light blue", size=.1) +
  #geom_density(alpha=0.1, fill="#FF6666") +
  # geom_vline(aes(xintercept=mean(phi, na.rm=T)),  color="black", linetype="dashed", size=1) +
  labs(y = NULL,
       x = "Fixed to Variable Cost Ratio") + 
       #x = expression(paste("Fixed to Variable Cost Ratio (",phi,")"))) + 
  scale_y_continuous(labels = comma, limits=c(0,6000)) +
  scale_x_continuous(labels = comma, limits=c(0,6000)) +
  theme_bw() +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2))  

q3<-ggplot(data=mydata,aes(x=rho)) + 
  geom_histogram(binwidth=240,color="black", fill = "light blue", size=.1) +
  #geom_density(alpha=0.1, fill="#FF6666")  + # Overlay with transparent density plot
  # geom_vline(aes(xintercept=mean(rho, na.rm=T)),  color="black", linetype="dashed", size=1) +
  scale_y_continuous(labels = comma,limits=c(0,8100)) +
  scale_x_continuous(labels = comma, limits=c(0,6000)) +
  labs(y = "Observations",
       x = "Requirements") +
       #x = expression(paste("Average Requirements (",rho,")"))) +
  theme_bw() +
  theme(axis.ticks = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 8),
       # axis.text.y = element_blank(), 
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2))  

#cool but maybe not so useful yet...
#ggplot(df, aes(x=density, fill=dsspGapI)) + geom_density(alpha=.25)

q4<-ggplot(data=mydata,aes(x=gamma)) + 
  geom_histogram(binwidth=0.25, color="black", fill = "light blue", size=.1) +
  #geom_density(alpha=1, fill="#FF6666")  + 
  scale_y_continuous(labels = comma,limits=c(0,8100)) +
  #scale_x_continuous(labels = comma, limits=c(-4,4)) +
  labs(y = NULL,
       x = expression(paste(gamma," Characteristic"))) + 
       #x = expression(paste(log(gamma)," Characteristic"))) + 
  theme_bw() +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  labels = trans_format("log10", math_format(10^.x))) 


q4<-ggplot(data=mydata,aes(x=gamma)) + 
  geom_histogram(binwidth = 30) +
  #geom_density(alpha=1, fill="#FF6666")  + 
  #scale_y_continuous(labels = comma,limits=c(0,8100)) +
  #scale_x_continuous(labels = comma, limits=c(-4,4)) +
  labs(y = NULL,
       x = expression(paste(gamma," Characteristic"))) + 
       #x = expression(paste(log(gamma)," Characteristic"))) + 
  theme_bw() +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #labels = trans_format("log10", math_format(10^.x))) 

q5<-ggplot(data=mydata,aes(x=supplyNodes)) + geom_histogram(binwidth=2, color="black", fill = "light blue", size=.1) + 
  theme_bw() + 
   labs(y = "Observations",
       x = "Supply Nodes") +
  scale_y_continuous(labels = comma,limits=c(0,15500)) +
  theme(axis.ticks = element_blank(), 
       # axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

q6<-ggplot(data=mydata,aes(x=demandNodes)) + geom_histogram(binwidth=2, color="black", fill = "light blue", size=.1) + 
  theme_bw() +
    labs(y = NULL,
       x = "Demand Nodes") +
   scale_y_continuous(labels = comma,limits=c(0,15500)) +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 



#save off histograms for paper -------------------------#
pdf(file="histogramSPhiRhoGamma.pdf",width=8,height=8)  #
grid.arrange(q5, q6, q1, q2, q3, q4, ncol=2)                    #  <------------------------------------------------------------------< save pdf images
dev.off()                                               #
#save off histograms for paper -------------------------#



#statitics for network characteristics -------------------------------------------------------

quickStat(mydata$density)
quickStat(mydata$phi)
quickStat(mydata$rho)
quickStat(mydata$noV)
quickStat(mydata$gamma)
quickStat(mydata$logGamma)
quickStat(mydata$totSupply)
quickStat(mydata$supplyNodes)
quickStat(mydata$demandNodes)
quickStat(mydata$dsspGap)
quickStat(mydata$dsspIT)


cor(mydata[mydata$phiBin==1,"phi"],mydata[mydata$phiBin==1,"dsspGap"])
cor(mydata[mydata$phiBin==2,"phi"],mydata[mydata$phiBin==2,"dsspGap"])

sink(file="correlation data.txt")
cor.test(mydata$totSupply,mydata$dsspIT,method="kendall")
cor.test(mydata$density,mydata$dsspIT,method="kendall")
cor.test(mydata$phi,mydata$dsspIT,method="kendall")
cor.test(mydata$rho,mydata$dsspIT,method="kendall")
cor.test(mydata$gamma,mydata$dsspIT,method="kendall")
cor.test(mydata$logGamma,mydata$dsspIT,method="kendall")
cor.test(mydata$supplyNodes,mydata$dsspIT,method="kendall")
cor.test(mydata$demandNodes,mydata$dsspIT,method="kendall")

cor.test(mydata$totSupply,mydata$dsspGap,method="kendall")
cor.test(mydata$density,mydata$dsspGap,method="kendall")
cor.test(mydata$phi,mydata$dsspGap,method="kendall")
cor.test(mydata$rho,mydata$dsspGap,method="kendall")
cor.test(mydata$gamma,mydata$dsspGap,method="kendall")
cor.test(mydata$logGamma,mydata$dsspGap,method="kendall")
cor.test(mydata$supplyNodes,mydata$dsspGap,method="kendall")
cor.test(mydata$demandNodes,mydata$dsspGap,method="kendall")
cor.test(mydata$dsspIT,mydata$dsspGap,method="kendall")
sink()



#data exploration....
summary(mydata$supplyNodes)  #~30%
summary(mydata$demandNodes)  #~20%

ggplot(data=mydata,aes(x=supplyNodes)) + geom_histogram(binwidth=2, color="black", fill = "light blue", size=.1) + 
  theme_bw() +
theme(panel.border = element_blank(), 
      axis.text = element_text(size = 7),
      axis.title.x = element_text(vjust=0.01),
      axis.title.y = element_text(vjust=0.2),
      axis.ticks = element_blank())

ggplot(data=mydata,aes(x=demandNodes)) + geom_histogram(binwidth=2, color="black", fill = "light blue", size=.1) + 
  theme_bw() +
theme(panel.border = element_blank(), 
      axis.text = element_text(size = 7),
      axis.title.x = element_text(vjust=0.01),
      axis.title.y = element_text(vjust=0.2),
      axis.ticks = element_blank())


q5<-ggplot(data=mydata,aes(x=demandNodes, y=dsspGap)) + geom_boxplot(color="black", fill = "light blue",  size=0.25) + 
  theme_bw() + 
labs(y = "DSSP Optimality Gap (%)", x="Demand Nodes") +
  theme_bw() +
    scale_y_continuous(labels = comma,limits=c(0,8100)) +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

q6<-ggplot(data=mydata,aes(x=supplyNodes, y=dsspGap)) + geom_boxplot(color="black", fill = "light blue",  size=0.25) +
  theme_bw() + 
labs(y = "DSSP Optimality Gap (%)", x="Supply Nodes") +
    scale_y_continuous(labels = comma,limits=c(0,8100)) +
 theme(axis.ticks = element_blank(), 
        #axis.text.y = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2)) 

mydata$demandNodes<-factor(mydata$demandNodes)
mydata$supplyNodes<-factor(mydata$supplyNodes)



#evaluate DSSP(1) 

table(mydata$dsspGapI )
summary(mydata$dsspGap)

#histogram for paper ----------------------------------------------------------------------- dsspGap histrogram
q1<-qplot(data=subset(mydata,dsspGap>0),x=dsspGap, binwidth=2) + theme_bw() +    
  labs(y = "Observations", x = "DSSP Optimality Gap (%)") +     
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


#break up parameters into percentiles and 50-quantiles; take mean of dsspGap and dsspIT  and plot (with error bars)--------------------------------------

#rho and dsspGap------------------------------------------------------------------
x<-quantile(mydata$rho,  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
x["0%"]<- x["0%"]*.9
mydata$rhoBins <- cut(mydata$rho, breaks = x)

dt<-data.table(mydata[,c("seed","rho","rhoBins","dsspGap","dsspIT")])
dt2<-data.frame(dt[,list(obs=length(seed),meanRho=mean(rho),meanGap=mean(dsspGap),meanIT=mean(dsspIT),sdGap=sd(dsspGap),seGap=stderr(dsspGap),seIT=stderr(dsspIT)),by=rhoBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g1<-ggplot(data=dt2,aes(x=meanRho,y=meanGap)) + 
  geom_point(size=2,shape=18) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,12.75)) +
  labs(y = "Mean DSSP Optimality Gap (%)", x = "Average Requirements \n\n (a)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
      #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

#rho and dsspIT -------------------------------------------
limits <- aes(ymax = meanIT + seIT, ymin=meanIT - seIT)

g2<-ggplot(data=dt2,aes(x=meanRho,y=meanIT)) + 
  geom_point(size=2,shape=18) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,25)) +
  labs(y = "DSSP Iterations", x = "Average Requirements \n\n (b)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
      #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  + 
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")


#for paper! ----------------------------------------#
pdf(file="rhoPerformance.pdf",width=8,height=4)     #
grid.arrange(g1, g2, ncol=2)                        # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#


#phi
x<-quantile(mydata$phi,  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
x
x["0%"]<- x["0.0%"]*.9
mydata$phiBins <- cut(mydata$phi, breaks = x)

dt<-data.table(mydata[,c("seed","phi","phiBins","dsspGap","dsspIT")])
dt2<-data.frame(dt[,list(obs=length(seed),meanPhi=mean(phi),meanGap=mean(dsspGap),meanIT=mean(dsspIT),sdGap=sd(dsspGap),seGap=stderr(dsspGap),seIT=stderr(dsspIT)),by=phiBins])
names(dt2)

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g3<-ggplot(data=dt2,aes(x=meanPhi,y=meanGap)) + 
  geom_point(size=2,shape=18) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,12.75)) +
  labs(y = "Mean DSSP Optimality Gap (%)", x = "Fixed to Variable Cost Ratio\n\n (a)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  + 
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

limits <- aes(ymax = meanIT + seIT, ymin=meanIT - seIT)

g4<-ggplot(data=dt2,aes(x=meanPhi,y=meanIT)) + 
  geom_point(size=2,shape=18) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,25)) +
  labs(y = "DSSP Iterations", x = "Fixed to Variable Cost Ratio\n\n (b)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  + 
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

#for paper! ----------------------------------------#
pdf(file="phiPerformance.pdf",width=8,height=4)     #
grid.arrange(g3, g4, ncol=2)                        # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#




x<-quantile(mydata$gamma,  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
x["0%"]<- x["0%"]*.9
mydata$gammaBins <- cut(mydata$gamma, breaks = x)

dt<-data.table(mydata[,c("seed","logGamma","gamma","gammaBins","dsspGap","dsspIT","rho","phi")])
dt2<-data.frame(dt[,list(obs=length(seed),meanLogGamma=mean(logGamma),meanPhi=mean(phi),meanRho=mean(rho),meanGamma=mean(gamma),meanGap=mean(dsspGap),meanIT=mean(dsspIT),sdGap=sd(dsspGap),seGap=stderr(dsspGap),seIT=stderr(dsspIT)),by=gammaBins])
names(dt2)

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g5<-ggplot(data=dt2,aes(x=meanGamma,y=meanGap)) + 
  geom_point(size=2,shape=18) + 
   geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,12.75)) +
  labs(y = "Mean DSSP Optimality Gap (%)", x = expression(atop(paste(gamma, " Characteristic"), paste("(a)")))) +
  theme(panel.border = element_blank(), 
         axis.text = element_text(size = 8),
         axis.title.x = element_text(vjust=0.01),
         axis.title.y = element_text(vjust=0.25),
         axis.ticks = element_blank()) +
   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))) +
   stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

limits <- aes(ymax = meanIT + seIT, ymin=meanIT - seIT)

g6<-ggplot(data=dt2,aes(x=meanGamma,y=meanIT)) + 
  geom_errorbar(width=0.0,size=0.25,limits) +
  geom_point(size=2,shape=18) +  
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,25)) +
  labs(y = "DSSP Iterations", x = expression(atop(paste(gamma, " Characteristic"), paste("(b)")))) +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
   stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")
 

#for paper! ----------------------------------------#
pdf(file="gammaPerformance.pdf",width=8,height=4)   #
grid.arrange(g5, g6, ncol=2)                        # <------------------------------------------------------------------------< save pdf images
dev.off()                                           #
#for paper! ----------------------------------------#



x<-quantile(mydata$density,  probs = c(seq(0,1,by=.02)),labels=c(seq(0,1,by=.02)))
x
x["0.0%"]<- x["0.0%"]*.9
mydata$densityBins <- cut(mydata$density, breaks = x)

dt<-data.table(mydata[,c("seed","density","densityBins","dsspGap","dsspIT")])
dt2<-data.frame(dt[,list(obs=length(seed),meanDensity=mean(density),meanGap=mean(dsspGap),meanIT=mean(dsspIT),sdGap=sd(dsspGap),seGap=stderr(dsspGap),seIT=stderr(dsspIT)),by=densityBins])
names(dt2)

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g7<-ggplot(data=dt2,aes(x=meanDensity,y=meanGap)) + 
  geom_point(size=2,shape=18) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,12.5)) +
   labs(y = "Mean DSSP Optimality Gap (%)", x = "Density\n\n (a)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  + 
   stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")

limits <- aes(ymax = meanIT + seIT, ymin=meanIT - seIT)

g8<-ggplot(data=dt2,aes(x=meanDensity,y=meanIT)) + 
  geom_point(size=2,shape=18) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,25)) +
  labs(y = "DSSP Iterations", x = "Density\n\n (b)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  + 
   stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")


#for paper! ------------------------------------------#
pdf(file="densityPerformance.pdf",width=8,height=4)   #
grid.arrange(g7, g8, ncol=2)                          # <------------------------------------------------------------------------< save pdf images
dev.off()                                             #
#for paper! ------------------------------------------#




#examine the lowest values for gamma ----------------------

test<-subset(mydata,gamma<0.01)[,c("seed","gamma","logGamma","rho","phi","density","dsspGap","totSupply")]
test$testGapI <- (test$gamma>.001)

aggregate(test$rho, list(test$testGapI), FUN = mean)
dt<-data.table(test)
(dt[,list(obs=length(seed),meanGamma=mean(gamma),meanGap=mean(dsspGap),meanRho=mean(rho),meanPhi=mean(phi),meanSupply=mean(totSupply)),by=testGapI])








# Evaluate DSSP(psi) ---------------------------------------------------------------------------------------------------------------

#used for paper in a table.......
table(mydata$dsspGapI )
table(mydata$pBestGapI )
summary(mydata$dsspGap)
summary(mydata$pBestGap)

table((subset(mydata,dsspGapI==1)$pBestGapI )

summary(subset(mydata,dsspGapI==1)$dsspGap)
summary(subset(mydata,dsspGapI==1)$pBestGap)

summary(subset(mydata,pBestGapI==1)$dsspGap)
summary(subset(mydata,pBestGapI==1)$pBestGap)
      
#..........................................
      
      
#prepare the DSSP and p-DSSP overlay trend graphs : good stuff! ----------------------------------------------------------------------      
    
      

#NOW FOR gamma ----------------------------------- 
      
x<-quantile(mydata$gamma,  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
x["0%"]<- x["0%"]*.9
mydata$gammaBins <- cut(mydata$gamma, breaks = x)
      
dt<-data.table(mydata[,c("seed","logGamma","gamma","gammaBins","dsspGap","pBestGap","dsspIT","rho","phi")])
dt2<-data.frame(dt[,list(meanGamma=mean(gamma),
                         DSSP=mean(dsspGap),
                         pDSSP=mean(pBestGap),
                         sdGap=sd(dsspGap),
                         seGap=stderr(dsspGap),
                         sdpGap=sd(pBestGap),
                         sepGap=stderr(pBestGap)),by=gammaBins])
  
names(dt2)      
tempMelt<-melt(dt2[,c("gammaBins","meanGamma","seGap","sepGap","DSSP","pDSSP")], id=c("gammaBins","meanGamma","seGap","sepGap"))        

tempMelt$se<-tempMelt$seGap      
tempMelt[tempMelt$variable == "pDSSP","se"]<-tempMelt[tempMelt$variable == "pDSSP","sepGap"]
      
tempMelt$variable <- factor(tempMelt$variable)      
levels(tempMelt$variable)[levels(tempMelt$variable)=="pDSSP"] <- "Parameterized DSSP"      
      
limits1 <- aes(ymax = value + se, ymin=value - se)   
      
names(tempMelt)[names(tempMelt)=="variable"]<-"Method"
      
g1<-ggplot(data=tempMelt,aes(x=meanGamma,y=value,group=Method)) + 
    geom_errorbar(limits1, width=0.0) +
    geom_point(size=2, aes(shape=Method, fill=Method)) +
    scale_y_continuous(labels = comma, limits=c(0,12.75)) +
    scale_shape_manual(values=c(18,24)) +    
    scale_fill_manual(values=c("black","white")) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))  +
    stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")  + 
    theme_bw() + 
    labs(y = "Mean DSSP Optimality Gap (%)", x = expression(atop(paste(gamma, " Characteristic"), paste("(d)")))) +
    theme(panel.border = element_blank(), 
              axis.text = element_text(size = 8),
              axis.title.x = element_text(vjust=0.01),
              axis.title.y = element_text(vjust=0.25),
              axis.ticks = element_blank(),
              legend.key = element_blank(),
              legend.title=element_blank(),
              legend.text = element_text(size = 12),
              legend.position="bottom"
              ) 
      
      
 #NOW FOR PHI -----------------------------------     
      
x<-quantile(mydata$phi,  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
x["0%"]<- x["0%"]*.9
mydata$phiBins <- cut(mydata$phi, breaks = x)
      
dt<-data.table(mydata[,c("seed","phi","phiBins","dsspGap","pBestGap")])
dt2<-data.frame(dt[,list(meanPhi=mean(phi),
                               DSSP=mean(dsspGap),
                               pDSSP=mean(pBestGap),
                               sdGap=sd(dsspGap),
                               seGap=stderr(dsspGap),
                               sdpGap=sd(pBestGap),
                               sepGap=stderr(pBestGap)),by=phiBins])
            
tempMelt<-melt(dt2[,c("phiBins","meanPhi","seGap","sepGap","DSSP","pDSSP")], id=c("phiBins","meanPhi","seGap","sepGap"))        
tempMelt$se<-tempMelt$seGap      
tempMelt[tempMelt$variable == "pDSSP","se"]<-tempMelt[tempMelt$variable == "pDSSP","sepGap"]
      
tempMelt$variable <- factor(tempMelt$variable)      
levels(tempMelt$variable)[levels(tempMelt$variable)=="pDSSP"] <- "Parameterized DSSP"      
limits1 <- aes(ymax = value + se, ymin=value - se)         
names(tempMelt)[names(tempMelt)=="variable"]<-"Method"
      
g2<-ggplot(data=tempMelt,aes(x=meanPhi,y=value,group=Method)) + 
        geom_errorbar(limits1, width=0.0) +
        geom_point(size=2, aes(shape=Method, fill=Method)) +
        scale_shape_manual(values=c(18,24)) +    
        scale_y_continuous(labels = comma, limits=c(0,12.75)) +
        scale_fill_manual(values=c("black","white")) +
        stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")  + 
        theme_bw() + 
        labs(y = "Mean DSSP Optimality Gap (%)", x = "Fixed to Variable Cost Ratio \n\n (c)") +
        theme(panel.border = element_blank(), 
              axis.text = element_text(size = 8),
              axis.title.x = element_text(vjust=0.01),
              axis.title.y = element_text(vjust=0.25),
              axis.ticks = element_blank(),
              legend.key = element_blank()
        )       
      
      
      
      
#now for density -------------------------------     
      
x<-quantile(mydata$density,  probs = c(seq(0,1,by=.05)),labels=c(seq(0,1,by=.05)))
x["0.0%"]<- x["0.0%"]*.9
mydata$densityBins <- cut(mydata$density, breaks = x)
      
dt<-data.table(mydata[,c("seed","density","densityBins","dsspGap","pBestGap")])
dt2<-data.frame(dt[,list(meanDensity=mean(density),
                               DSSP=mean(dsspGap),
                               pDSSP=mean(pBestGap),
                               sdGap=sd(dsspGap),
                               seGap=stderr(dsspGap),
                               sdpGap=sd(pBestGap),
                               sepGap=stderr(pBestGap)),by=densityBins])
      
tempMelt<-melt(dt2[,c("densityBins","meanDensity","seGap","sepGap","DSSP","pDSSP")], id=c("densityBins","meanDensity","seGap","sepGap"))        
tempMelt$se<-tempMelt$seGap      
tempMelt[tempMelt$variable == "pDSSP","se"]<-tempMelt[tempMelt$variable == "pDSSP","sepGap"]
tempMelt$variable <- factor(tempMelt$variable)      
levels(tempMelt$variable)[levels(tempMelt$variable)=="pDSSP"] <- "Parameterized DSSP"      
limits1 <- aes(ymax = value + se, ymin=value - se)         
names(tempMelt)[names(tempMelt)=="variable"]<-"Method"
      
g3<-ggplot(data=tempMelt,aes(x=meanDensity,y=value,group=Method)) + 
        geom_errorbar(limits1, width=0.0) +
        geom_point(size=2, aes(shape=Method, fill=Method)) +
        scale_shape_manual(values=c(18,24)) +    
        scale_y_continuous(labels = comma, limits=c(0,12.75)) +
        scale_fill_manual(values=c("black","white")) +
        stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")  + 
        theme_bw() + 
        labs(y = "Mean DSSP Optimality Gap (%)", x = "Density \n\n (b)") +
        theme(panel.border = element_blank(), 
              axis.text = element_text(size = 8),
              axis.title.x = element_text(vjust=0.01),
              axis.title.y = element_text(vjust=0.25),
              axis.ticks = element_blank(),
              legend.key = element_blank()
        )         
      
      
      
#now for Rho -------------------
      
x<-quantile(mydata$rho,  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
x["0%"]<- x["0%"]*.9
mydata$rhoBins <- cut(mydata$rho, breaks = x)
      
dt<-data.table(mydata[,c("seed","rho","rhoBins","dsspGap","pBestGap")])
dt2<-data.frame(dt[,list(meanRho=mean(rho),
                               DSSP=mean(dsspGap),
                               pDSSP=mean(pBestGap),
                               sdGap=sd(dsspGap),
                               seGap=stderr(dsspGap),
                               sdpGap=sd(pBestGap),
                               sepGap=stderr(pBestGap)),by=rhoBins])
      
tempMelt<-melt(dt2[,c("rhoBins","meanRho","seGap","sepGap","DSSP","pDSSP")], id=c("rhoBins","meanRho","seGap","sepGap"))        
tempMelt$se<-tempMelt$seGap      
tempMelt[tempMelt$variable == "pDSSP","se"]<-tempMelt[tempMelt$variable == "pDSSP","sepGap"]
tempMelt$variable <- factor(tempMelt$variable)      
levels(tempMelt$variable)[levels(tempMelt$variable)=="pDSSP"] <- "Parameterized DSSP"      
limits1 <- aes(ymax = value + se, ymin=value - se)         
names(tempMelt)[names(tempMelt)=="variable"]<-"Method"
      
g4<-ggplot(data=tempMelt,aes(x=meanRho,y=value,group=Method)) + 
        geom_errorbar(limits1, width=0.0) +
        geom_point(size=2, aes(shape=Method, fill=Method)) +
        scale_shape_manual(values=c(18,24)) +  
        scale_y_continuous(labels = comma, limits=c(0,12.75)) +
        scale_fill_manual(values=c("black","white")) +
        stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")  + 
        theme_bw() + 
        labs(y = "Mean DSSP Optimality Gap (%)", x = "Requirements \n\n (a)") +
        theme(panel.border = element_blank(), 
              axis.text = element_text(size = 8),
              axis.title.x = element_text(vjust=0.01),
              axis.title.y = element_text(vjust=0.25),
              axis.ticks = element_blank(),
              legend.key = element_blank()
        )       
      
      
      
            
mylegend<-g_legend(g1)
      
#for paper! --------------------------------------------------------#
pdf("EvalPSi.pdf",width=8,height=8)                                 #
grid.arrange(arrangeGrob(g4 + theme(legend.position="none"),        #
                               g3 + theme(legend.position="none"),  #
                               g2 + theme(legend.position="none"),  # <---------------------------------------------------------------< save pdf images
                               g1 + theme(legend.position="none"),  #
                               nrow=2),                             #    
                         mylegend, nrow=2,heights=c(10, 1))         # 
dev.off()                                                           #
#for paper! --------------------------------------------------------#      
      
      
      
      
cor.test(mydata$totSupply,mydata$pBestGap)$estimate
cor.test(mydata$density,mydata$pBestGap)$estimate
cor.test(mydata$phi,mydata$pBestGap)$estimate
cor.test(mydata$rho,mydata$pBestGap)$estimate
cor.test(mydata$gamma,mydata$pBestGap)$estimate
cor.test(mydata$logGamma,mydata$pBestGap)$estimate
cor.test(mydata$supplyNodes,mydata$pBestGap)$estimate
cor.test(mydata$demandNodes,mydata$pBestGap)$estimate
   

skinny<-mydata[,c("dsspGap","pBestGap","totSupply","phi","rho","logGamma","supplyNodes","demandNodes","gamma","density")]      
      
group2 <- skinny[sample(1:nrow(skinny), 50, replace=TRUE),]

      
R = 1000                                        # the number of replicates
pBestCor <- numeric(R)
dsspCor <- numeric(R)      
corr.values <- data.frame(pBestCor,dsspCor)                        # storage for the results
for (i in 1:R) {
group1 <- skinny[sample(1:nrow(skinny), 1000, replace=TRUE),]
df<-mydata[sample(nrow(mydata),100),]
corr.values$pBestCor[i] <- cor.test(group1$phi,group1$pBestGap)$estimate
group2 <- skinny[sample(1:nrow(skinny), 1000, replace=TRUE),]
corr.values$dsspCor[i] <- cor.test(group1$phi,group1$dsspGap)$estimate
}
t.test(group1,group2)$statistic

      
t.test(mydata$dsspGap,mydata$pBestGap)
      
      

  
bootStrapCorr("totSupply")
bootStrapCorr("supplyNodes")
bootStrapCorr("demandNodes")
bootStrapCorr("density")
bootStrapCorr("phi")
bootStrapCorr("rho")
bootStrapCorr("gamma")
bootStrapCorr("logGamma")



      
      
#histograms for psi   -------------- good stuff
      

      
dt<-data.table(mydata[,c("seed","pBest","pBest2")])
tempMelt<-melt(dt, id=c("seed"))        
tempMelt$variable<-factor(tempMelt$variable)

levels(tempMelt$variable)[levels(tempMelt$variable)=="pBest"]<-"Minimum"
levels(tempMelt$variable)[levels(tempMelt$variable)=="pBest2"]<-"Maximum"
      
      

      
g1<-ggplot(tempMelt, aes(x=value)) + facet_wrap(~variable) +  
        geom_histogram(binwidth=0.1, color="black", fill = "light blue", size=.1) +   
        theme_bw() +
        labs(y = "Observations", x = "DSSP Parameter Value") +
        theme(panel.border = element_blank(), 
              axis.text = element_text(size = 8),
              axis.title.x = element_text(vjust=0.01),
              axis.title.y = element_text(vjust=0.25),
              axis.ticks = element_blank(),
              legend.key = element_blank())      
    

pdf("psiHistograms.pdf",width=8,height=4)
print (g1)
dev.off()     
      
      
      
      
summary(mydata$pBest)
summary(mydata$pBest2)
      
      
table(mydata$pBestIL)
table(mydata$pBest2IL)
    
     
      
table(mydata$pBestBetterI)      
     
table(subset(mydata,dsspGapI==1)$pBestBetterI)    
table(subset(mydata,pBestIL==1)$pBestBetterI)      
table(subset(mydata,pBest2IL==1)$pBestBetterI)        
      
subset(mydata,pBest2IL==1 & seed < 100)[,c("seed","dsspGap","pBestGap","pBest")] 
      
     
      
      
summary(subset(mydata,pBestIL==1 & pBestBetterI ==1)$dsspGap)      
summary(subset(mydata,pBestIL==1 & pBestBetterI ==1)$pBestGap)        
      
summary(subset(mydata,pBest2IL==1 & pBestBetterI ==1)$dsspGap)
summary(subset(mydata,pBest2IL==1 & pBestBetterI ==1)$pBestGap)    
      
      
      
      

mydata$density50 <- (mydata$density > 50) 
mydata$phi1000 <- (mydata$phi > 50) 
mydata$gamma0.05 <- (mydata$gamma > 0.050)       

mydata$supplyNodesV <- as.numeric(levels(mydata$supplyNodes))[mydata$supplyNodes]/12
mydata$demandNodesV <- as.numeric(levels(mydata$demandNodes))[mydata$demandNodes]/12

      
      
mydata$pbestGapT<-mydata$pBestGap #log(mydata$pBestGap+1)      
mydata$dsspGapT<-mydata$dsspGap #log(mydata$dsspGap+1)     
mydata$ImprovementT <- mydata$dsspGapT - mydata$pbestGapT

      
trimData<-subset(mydata,dsspGapI==1)      
summary(trimData)      
      
lfit<-lm(data=trimData, ImprovementT ~ logGamma * gamma0.05 + phi * phi1000 + rho + density + supplyNodesV + demandNodesV + totSupply)
summary(lfit)
plot(lfit)
    


lfitD<-lm(data=trimData, dsspGapT ~ logGamma * gamma0.05 + phi * phi1000 + rho + density + supplyNodesV + demandNodesV + totSupply)
summary(lfitD)     
 
lfitP<-lm(data=trimData, pbestGapT ~ logGamma * gamma0.05 + phi * phi1000 + rho + density + supplyNodesV + demandNodesV + totSupply)
summary(lfitP)        
      
      
df<-data.frame(predD=lfitD$fitted.values, residD=lfitD$residuals, actualD=lfitD$fitted.values+lfitD$residuals,
               predP=lfitP$fitted.values, residP=lfitP$residuals, actualP=lfitP$fitted.values+lfitP$residuals,
               difFit=lfitD$fitted.values-lfitP$fitted.values,
               difAct=lfitD$fitted.values+lfitD$residuals - lfitP$fitted.values -lfitP$residuals)        


df$pPrank<-rank(df$predP, ties.method= "random")
df$pDrank<-rank(df$predD, ties.method= "random")
df$aPrank<-rank(df$actualP, ties.method= "random")
df$aDrank<-rank(df$actualD, ties.method= "random")

df$actualPDeciles <- cut(df$aPrank, 10, include.lowest = TRUE,labels=c(10:1))       
df$actualDDeciles <- cut(df$aDrank, 10, include.lowest = TRUE,labels=c(10:1))       
df$predPDeciles <- cut(df$pPrank, 10, include.lowest = TRUE,labels=c(10:1))       
df$predDDeciles <- cut(df$pDrank, 10, include.lowest = TRUE,labels=c(10:1))       

table(df$predDDeciles,df$actualDDeciles)

      
      
dt<-data.table(df)
dt2<-data.frame(dt[,list(count=length(seed),by=c(df$predDDeciles,df$actualDDeciles))])
    
      
      x<-quantile(df$predP,  probs = c(seq(0,1,by=.1)),labels=c(seq(0,1,by=.1)))
x
df$predPDeciles <- cut(df$predP, breaks = x, include.lowest = TRUE) 
x<-quantile(df$actualP,  probs = c(seq(0,1,by=.1)),labels=c(seq(0,1,by=.1)))
x
df$actualPDeciles <- cut(df$actualP, breaks = x, include.lowest = TRUE) 
      
      
      
qplot(data=df, x=predPDeciles), y=difFit, alpha=I(1/10))      

qplot(data=df, x=predD), y=difFit, alpha=I(1/10))      






#map parameter space

MAPpsi <- read.csv("~/Research and Grants/Research/Time Varying Cost/DSSP search results/MAPpsi.txt")
names(MAPpsi)

table(MAPpsi$seed)
#MAPpsi$gap<-MAPpsi$gap*100

MAPpsi$seed.f<-factor(MAPpsi$seed)
levels(MAPpsi$seed.f)<-c(paste("Instance",c("11" ,"12" ,"16", "20" ,"29" ,"49" ,"60")))

g1<-ggplot(data=subset(MAPpsi, seed== 11 | seed == 20 | seed == 60),aes(x=pVal,y=gap,group=seed)) + 
    facet_wrap(~ seed.f, ncol=1) + 
    geom_line(size=0.65) + 
    theme_bw() +
    labs(y = "Optimality Gap (%)", x = "DSSP Parameter Value") +
        theme(panel.border = element_blank(), 
              axis.text = element_text(size = 8),
              axis.title.x = element_text(vjust=0.01),
              axis.title.y = element_text(vjust=0.25),
              axis.ticks = element_blank(),
              legend.key = element_blank())    

pdf("psiSpace.pdf",width=5,height=8)
print(g1)
dev.off()

MAPpsi1<-subset(MAPpsi,pVal==1)

MAPmerge<-merge(MAPpsi,MAPpsi1[,c("seed","pObj")],by="seed")
MAPmerge$BetterThanPsi1 <- (MAPmerge$pObj.x <= MAPmerge$pObj.y)
by(MAPmerge$BetterThanPsi1, MAPmerge$seed, table)
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
tempMelt<-melt(mydata[,c("seed","dsspGap","pBestGap","dsspGapI","pBestGapI")], id=c("seed","dsspGapI","pBestGapI"))
      
   

min(mydata[mydata$pBestGapI,]$pBestGap)

ggplot(df, aes(x=density, fill=dsspGapI)) + geom_density(alpha=.25)


      
      
#histogram for paper ----------------------------------------------------------------------- dsspGap histrogram
q1<-qplot(data=subset(mydata,pBestGap>0),x=pBestGap,binwidth=2) + theme_bw() +    
  labs(y = "Observations", x = "Parameterized DSSP Optimality Gap (%)") +     
  geom_histogram(binwidth=2, color="black", fill = "light blue", size=.1) +   
  scale_y_continuous(labels = comma, limits=c(0,3500)) +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank()
  )  

q2<-qplot(data=subset(mydata,pBestGap>0),x=dsspGap,binwidth=2) + theme_bw() +    
  labs(y = "Observations", x = "DSSP Optimality Gap (%)") +     
  geom_histogram(binwidth=2, color="black", fill = "light blue", size=.1) +   
  scale_y_continuous(labels = comma, limits=c(0,3500)) +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank()
  ) 

q3<-qplot(data=subset(mydata,dsspGap>0),x=pBestGap,binwidth=2) + theme_bw() +    
  labs(y = "Observations", x = "Parameterized DSSP Optimality Gap (%)") +     
  geom_histogram(binwidth=2, color="black", fill = "light blue", size=.1) +   
  scale_y_continuous(labels = comma, limits=c(0,3500)) +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank()
  )  

q4<-qplot(data=subset(mydata,dsspGap>0),x=dsspGap,binwidth=2) + theme_bw() +    
  labs(y = "Observations", x = "DSSP Optimality Gap (%)") +     
  geom_histogram(binwidth=2, color="black", fill = "light blue", size=.1) +   
  scale_y_continuous(labels = comma, limits=c(0,3500)) +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank()
  ) 


grid.arrange(q1, q2, q3, q4, ncol=2)




qplot(data=tempMelt, x=value, fill=variable, facets=pBestGapI~dsspGapI)


ggplot(tempMelt, aes(x=value, fill=variable)) +   facet_wrap(dsspGapI~pBestGapI) +
  geom_density(alpha=.75) +
  scale_fill_brewer(palette="Blues") 



qplot(data=mydata,x=dsspGap,y=pBestGap,geom="point",alpha=I(1/10))

+scale_fill_brewer(palette="Blues")

mydata$dsspGapI<-factor(mydata$dsspGapI)

ggpairs(mydata, columns=c("dsspGap","pBestGap","dsspGapI"), 
        diag=list(continuous="density",   discrete="bar"), axisLabels="show")

#rho and dsspGap------------------------------------------------------------------
x<-quantile(mydata$rho,  probs = c(seq(0,1,by=.01)),labels=c(seq(0,1,by=.01)))
x["0%"]<- x["0%"]*.9
mydata$rhoBins <- cut(mydata$rho, breaks = x)

dt<-data.table(mydata[,c("seed","rho","rhoBins","dsspGap","dsspIT")])
dt2<-data.frame(dt[,list(obs=length(seed),meanRho=mean(rho),meanGap=mean(dsspGap),meanIT=mean(dsspIT),sdGap=sd(dsspGap),seGap=stderr(dsspGap),seIT=stderr(dsspIT)),by=rhoBins])

limits <- aes(ymax = meanGap + seGap, ymin=meanGap - seGap)

g1<-ggplot(data=dt2,aes(x=meanRho,y=meanGap)) + 
  geom_point(size=2,shape=18) + 
  geom_errorbar(limits, width=0.0) +
  theme_bw() +
  scale_y_continuous(labels = comma, limits=c(0,12.5)) +
  labs(y = "Mean DSSP Optimality Gap (%)", x = "Average Requirements \n\n (a)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  +
  stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")















ggplot(data=test,aes(x=gamma,y=rho,color=testGapI)) + 
  geom_point(size=4) + 
  theme_bw() +
  #scale_y_continuous(labels = comma, limits=c(0,12.5)) +
  #labs(y = "Mean DSSP Optimality Gap (%)", x = "Density\n\n (a)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 8),
        #  axis.title = element_text(size = 8),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.ticks = element_blank())  + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
stat_smooth(method = "loess",size=0.5,se=FALSE,color="black")




qplot(data=mydata,x=densityBins,y=dsspIT,geom="boxplot")



mydata$demandNodesF <- factor(mydata$demandNodes)
mydata$supplyNodesF <- factor(mydata$supplyNodes)


fit <- lm(data=mydata,dsspGap~density+rho+phi+logGamma+demandNodesF+supplyNodesF)
summary(fit)

ggplot(data=mydata,aes(y=dsspGap, x=demandNodesF)) + facet_wrap(~supplyNodesF) +
  geom_boxplot(outlier.size=.5, color="black", fill = "light blue",  size=0.25) +
  theme_bw() +
  labs(y = "DSSP Iterations", x = "Total Network Supply \n\n (a)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank())  






aggregate(mydata$dsspGap, list(mydata$rhoBins), FUN = median)
aggregate(mydata$dsspGap, list(mydata$densityBins), FUN = median)
aggregate(mydata$dsspGap, list(mydata$densityBins), FUN = quantile)

quantile(mydata[mydata$density,"dsspGap"], seq(0,1,by=.1))

aggregate(mydata$dsspGap, list(mydata$rhoBins), FUN = median)
aggregate(mydata$dsspGap, list(mydata$rhoBins), FUN = median)











g2<-ggplot(data=mydata,aes(x=logGamma, y=dsspGap)) +
  geom_point( size=0.75, alpha=I(1/30)) +
  theme_bw() + coord_cartesian(ylim = c(-3,45)) +
  stat_smooth(fill="blue", color="darkblue", size=1,alpha=0.3) +
  labs(y = "DSSP Optimality Gap (%)", x=expression(log(gamma))) +
  theme(panel.border = element_blank(), 
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.text = element_text(size = 7),
        axis.ticks = element_blank()) 
       

pdf(file="gammaTrend.pdf",width=8,height=3.5)
print(g2)
dev.off()


pdf(file="gamma.pdf",width=8,height=7)
grid.arrange(g1, g2, ncol=1)
dev.off()




%dsspGap ~ logGamma + phi^3? + phi^2 + phi + totSupply + rho + density + supplyNodes + demandNodes + phiExtreme

mydata$phi2 <- mydata$phi*mydata$phi
mydata$phi3 <- mydata$phi*mydata$phi2


lfit<-lm(data=mydata, dsspIT ~  logGamma + gammaBins + phi + density + supplyNodes + demandNodes)
summary(lfit)
plot(lfit)


lfit<-aov(data=mydata, dsspGap ~  gammaBins)
plot(TukeyHSD(lfit, "gammaBins",ordered = TRUE))

lfit<-aov(data=mydata, dsspGap ~  rhoBins)
plot(TukeyHSD(lfit, "rhoBins",ordered = TRUE))

lfit<-lm(data=subset(mydata,phi<1000), dsspGap ~ phi2 + phi + logGamma + density + logGamma * gammaBins +  supplyNodes + demandNodes)
summary(lfit)

plot(mydata$phi,fitted(lfit))


lfit<-lm(data=subset(mydata,phi<1000), dsspGap ~ logGamma * gammaBins + totSupply + phi + rho + density + supplyNodes + demandNodes)
summary(lfit)

step <- stepAIC(lfit, direction="both")
step$anova # display results

lfit<-lm(data=subset(mydata,dsspGap>0), dsspGap ~ logGamma * gammaBins + totSupply +  rho + density + supplyNodes + demandNodes)
summary(lfit)
str(lfit)

df<-data.frame(actual=fitted(lfit)+residuals(lfit),fitted=fitted(lfit),resid=residuals(lfit))
qplot(data=df, x=fitted, y=resid, alpha=I(1/20))

dev.off()


mydata$dsspGapI

lrfitM<-glm(data=mydata,dsspGapI~logGamma * gammaBins + totSupply +  rho + density + supplyNodes + demandNodes)
summary(lrfitM)


%dsspIT?

lfit<-lm(data=mydata, dsspIT ~ logGamma + density + phi + rho  + supplyNodes + demandNodes + gammaBins)
summary(lfit)


qplot(data=mydata,x=rho,y=dsspGap,alpha=I(1/10))
qplot(data=subset(mydata[mydata$dsspGapI,],phi<300),x=phi,y=dsspGap,alpha=I(1/5),geom=(c("point","smooth")))

qplot(data=mydata,x=density,y=dsspGap,alpha=I(1/10)) + geom_smooth(method="lm")
qplot(data=mydata, x=densityBins, y=dsspGap, geom="boxplot")
qplot(data=mydata, x=rhoBins, y=dsspGap, geom="boxplot")

qplot(x=bins10,  y=dsspGap, data=subset(mydata,phi<100), geom="boxplot", facets = densityBins ~ ., outline=FALSE)

qplot(x=bins,  y=dsspGap, data=subset(mydata,phi<1000), geom="boxplot", facets = densityBins ~ .)


qplot(data=mydata,x=gamma,y=dsspGap,alpha=I(1/5)) + theme_bw() + 
  scale_x_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))

lfit<-lm(data=mydata,dsspGap~factor(rhoBins))
summary(lfit)

x<-rnorm(500)
df<-data.frame(x)



q<-quantile(x,  probs = c(seq(0,1,by=.25)))
q
1.5*IQR(x)

a_plot.labels <-data.frame(
  yval = c(-0.75,-0.0918,0.64,-0.75-2.037,0.64+2.037),
  xval = c(1.2,1.2,1.2,1.1,1.1),
  label = c("Q1","median","Q3","Q1 - 1.5 x IQR","Q3 + 1.5 x IQR")
)

a_plot <- ggplot(data=df, aes(factor(0),x)) + geom_boxplot(color="black", fill = "light blue", outlier.shape = NA, width=0.35)+ coord_flip(ylim=c(-3.5,3.5)) + 
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid=element_blank()) + 
  geom_text(data = a_plot.labels, aes(x = xval, y = yval, label = label),size=4)

a_plot


#A viewport taking up a fraction of the plot area
vp <- viewport(width = 0.35, height = 0.35, x = 0.8, y = 0.8)

pdf(file="test.pdf",height=6,width=4)
#Just draw the plot twice
print(b)
print(a_plot, vp = vp)
dev.off()




p = qplot(1:10, 1:10, log="x")
g = qplot(0.1 , 0.1)
vp1 <- viewport(width = 0.5, 
                height = 0.3, 
                x = 0.4, 
                y = 0.7)
vp2 <- viewport(width = 0.3, 
                height = 0.3, 
                x = 0.8, 
                y = 0.3)
#png("text.png")
print(p)
print(g, vp = vp1)
print(a_plot, vp = vp2)









# create boxplot that includes outliers
qplot(x=bins50,  y=dsspGap, data=subset(mydata,phi<1000), geom="boxplot")

# create boxplot that includes outliers
p0=qplot(x=bins50,  y=dsspGap, data=subset(mydata,phi<1000), geom="boxplot",outlier.shape = NA)

# compute lower and upper whiskers
ylim1 = boxplot.stats(subset(mydata,phi<1000)$dsspGap)$stats[c(1, 5)]
ylim1

# scale y limits based on ylim1
p1 = p0 + coord_cartesian(ylim = 1.5*ylim1)
p1

myPlots <- function(DS, x, y,ylim1) {
  aes <- eval(substitute(aes(x, y),
                         list(x = substitute(x), y = substitute(y))))
  p <- ggplot(DS, aes) + geom_boxplot(outlier.shape = NA) 
  p1 = p + coord_cartesian(ylim = ylim1)
  
}

p<-myPlots(subset(mydata,phi<1000),bins50,dsspGap,c(0,40))
p
p<-myPlots(subset(mydata,phi<1000),gammaBins,pBest,c(0,1.25))
p
p<-myPlots(subset(mydata,phi<1000),rhoBins,dsspGap,c(0,40))
p




df<-mydata[sample(nrow(mydata),100),]

ggpairs(df, columns=c("rhoBins", "densityBins", "gammaBins", "dsspGap"), 
        diag=list(continuous="density",   discrete="bar"), axisLabels="show")

ggpairs(df, columns=c("pBestLog", "phi"), 
        diag=list(continuous="density"), axisLabels="show")

qqplot(data=df,x=pBestlog,y=phi)

boxplot.stats(mydata["dsspGap"])$stats[c(1, 5)]
ylim1 = boxplot.stats(mydata$dsspGap)

# create boxplot that includes outliers
p0 = qplot(x=bins10,  y=dsspGap, data=subset(mydata,phi<100), geom="boxplot", facets = densityBins ~ .,outlier.shape = NA)
p0

# compute lower and upper whiskers
ylim1 = boxplot.stats(df$dsspGap)$stats[c(1, 5)]
ylim1 = c(0,13)

# scale y limits based on ylim1
p1 = p0 + coord_cartesian(ylim = ylim1)
p1



qplot(data=df,x=pBestLog,y=phi)







qplot(data=mydata,x=rho,y=noV)
qplot(data=mydata,x=density,y=pBest)


qplot(data=mydata, x=densityBins, y=pBest, geom="boxplot",outlier.shape = NA)
qplot(data=mydata, x=densityBins, y=dsspGap, geom="boxplot")

qplot(x=bins10,  y=dsspGap, data=subset(mydata,phi<100), geom="boxplot", facets = densityBins ~ ., outline=FALSE)

qplot(x=bins,  y=dsspGap, data=subset(mydata,phi<1000), geom="boxplot", facets = densityBins ~ .)






qplot(data=mydata,x=phi,y=pBest)

qplot(data=mydata,x=gammaBins,y=pBest,geom="boxplot")

qplot(data=mydata,x=rho,y=dsspGap)
qplot(data=mydata,x=phi,y=dsspGap)
qplot(data=mydata,x=gamma,y=dsspGap) +
  scale_x_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))

qplot(data=mydata,x=gammaBins,y=dsspGap,geom="boxplot")

qplot(x=bins,  y=dsspGap, data=subset(mydata,phi<1000), geom="boxplot", facets = rhoBins ~ .)
qplot(x=bins,  y=dsspGap, data=subset(mydata,phi<1000), geom="boxplot", facets = gammaBins ~ .)

qplot(data=mydata, x=rhoBins)


# create boxplot that includes outliers
p0 = qplot(x=bins10,  y=dsspGap, data=subset(mydata,phi<100), geom="boxplot", facets = densityBins ~ .,outlier.shape = NA)
p0

# compute lower and upper whiskers
ylim1 = boxplot.stats(df$dsspGap)$stats[c(1, 5)]
ylim1 = c(0,13)

# scale y limits based on ylim1
p1 = p0 + coord_cartesian(ylim = ylim1)
p1


#learning <- as density increases; dsspGap increases

df<-mydata[sample(nrow(mydata),2000),]

ggpairs(df, columns=c("gammaBins", "densityBins", "dsspGap", "dsspIT"), 
        diag=list(continuous="density",   discrete="bar"), axisLabels="show")



ggpairs(df, columns=c("supplyNodesF", "demandNodesF", "dsspGap", "dsspIT"), 
        diag=list(continuous="density",   discrete="bar"), axisLabels="show")


ggpairs(
  tips[,1:4],
  upper = list(continuous = "density", combo = "facethist"),
  lower = list(continuous = "points", combo = "box")
)


#some old stuff not using now, but dont want to erase it yet----------------------------------------------

#mydata$bins10 <- cut(mydata$phi, breaks = c(0, seq(10, 5400, by = 10)), labels = 0:539)
#mydata$phiBinsQ <- cut(mydata$phi, breaks = c(0, 500, 1000, 1500, 2000), labels = 0:3)
#mydata$phiBins <- cut(mydata$phi, breaks = c(0, seq(100, 2000, by = 100)), labels = 0:19)
#mydata$bins50 <- cut(mydata$phi, breaks = c(0, seq(50, 2000, by = 50)), labels = 0:39)
#mydata$PhiExtreme <- cut(mydata$phi, breaks = c(0,30, 1950,2000),labels=0:2)

#binomial distributions; not needed for paper
#qplot(data=mydata,supplyNodes,stat="bin") + scale_x_discrete("Supply Nodes") + ylab("Observations")
#qplot(data=mydata,demandNodes,stat="bin") + scale_x_discrete("Supply Nodes") + ylab("Observations")


#pdf(file="histograms1.pdf",width=8,height=4)
#grid.newpage()
#pushViewport(viewport(layout = grid.layout(1, 2)))
#print(q1, vp = vplayout(1, 1))
#print(q2, vp = vplayout(1, 2))
#dev.off()
#dev.off()


#for paper! -----------------------------------------------------------------<
gS<-ggplot(data=mydata,aes(y=dsspGap, x=totSupply)) + 
  geom_point(alpha=I(1/40),size=0.7) +  
  stat_smooth(method="lm",formula=y~x,fill="light blue", color="darkblue", size=1,alpha=0.5) + 
  theme_bw() +
  labs(y = "DSSP Optimality Gap (%)", x = "Total Network Supply \n\n (a)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank())  


#for paper! -----------------------------------------------------------------<
gD<-ggplot(data=mydata,aes(y=dsspGap, x=densityBins)) + 
  geom_boxplot(outlier.size=.5, color="black", fill = "light blue",  size=0.25) +
  theme_bw() +
  labs(y = "DSSP Optimality Gap (%)", x = "Graph Density (%) \n\n (b)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank()
  )  

#for paper! -----------------------------------------------------------------<
gR<-ggplot(data=mydata,aes(y=dsspGap, x=rhoBins)) + 
  geom_boxplot(outlier.size=.5, color="black", fill = "light blue",  size=0.25) +
  theme_bw() +
  labs(y = "DSSP Optimality Gap (%)", x = "Average Requirement \n\n (c)") +
  theme(panel.border = element_blank(), 
        axis.text = element_text(size = 7),
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.2),
        axis.ticks = element_blank()
  )  


gP<-ggplot(data=mydata,aes(x=phi,y=dsspGap))  + 
  stat_smooth(method="lm",formula=y~poly(x,3),fill="blue", color="darkblue", size=1,alpha=0.5) +
  theme_bw() +
  labs(y = "DSSP Optimality Gap (%)", x="Fixed to Variable Cost Ratio Estimate \n\n (d)") +
  theme(panel.border = element_blank(), 
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.text = element_text(size = 7),
        axis.ticks = element_blank()) +
  stat_summary(fun.y = "mean", geom = "point", shape= 17, size= 1.2, color="darkblue",fill= "blue",alpha=I(1/30))


#pdf(file="relationships.pdf",width=8,height=8)
grid.arrange(gS, gD, gR, gP, ncol=2)
#dev.off()


pdf(file="phiTrend.pdf",width=6,height=5)
ggplot(data=mydata,aes(x=phi,y=dsspGap))  + 
  stat_smooth(method="lm",formula=y~poly(x,3),fill="blue", color="darkblue", size=1,alpha=0.5) +
  theme_bw() +
  labs(y = "DSSP Optimality Gap (%)", x=expression(phi)) +
  theme(panel.border = element_blank(), 
        axis.title = element_text(size=9) , 
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6),
        axis.ticks = element_blank()) +
  stat_summary(fun.y = "mean", geom = "point", shape= 20, size= 1.5, fill= "black",alpha=I(1/8))
dev.off()


pdf(file="phiBins.pdf",width=7,height=4)
ggplot(data=mydata,aes(x=phiBins, y=dsspGap)) + 
  geom_boxplot(outlier.size=.5, color="black", fill = "light blue",  size=0.25) + 
  theme_bw() + coord_cartesian(ylim = c(0,70)) +
  labs(y = "DSSP Optimality Gap (%)", x=expression(paste(gamma, " deciles"))) +
  theme(panel.border = element_blank(), 
        axis.title = element_text(size=9) , 
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6),
        axis.ticks = element_blank()) +
  stat_summary(fun.y = "mean", geom = "point", shape= 18, size= 1.5, fill= "black")
dev.off()

#---for paper ---------------------------------------------------------------------<
x<-quantile(mydata$gamma,  probs = c(seq(0,1,by=.1)),labels=c(seq(0,1,by=.1)))
x
x["0%"]<- 0

#good stuff
mydata$gammaBins <- cut(mydata$gamma, breaks = x)
mydata$gammaBins <- cut(mydata$gamma, breaks = x, labels=paste(labels(x)[seq(1,length(x)-1)],"-",labels(x)[seq(2,length(x))]))

my_labs <- ggplot2:::interleave(seq(1,100,by=2), "")
mydata$gammaBins <- cut(mydata$gamma, breaks = x, labels=seq(1,length(x)-1))


levels(mydata$gammaBins)[10]<-"(3.77,14580]"

g1<-ggplot(data=mydata,aes(x=gammaBins, y=dsspGap)) + geom_boxplot(outlier.shape = NA, color="black", fill = "light blue",  size=0.25) + 
  theme_bw() + coord_cartesian(ylim = c(0,45)) + facet_wrap(~ rhoBins)
labs(y = "DSSP Optimality Gap (%)", x=expression(paste(gamma, " deciles"))) +
  theme(panel.border = element_blank(), 
        axis.title.x = element_text(vjust=0.01),
        axis.title.y = element_text(vjust=0.25),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.ticks = element_blank())

qplot(data=mydata,x=phi,y=gamma, geom="smooth")


pdf(file="gammaDeciles.pdf",width=8,height=3.4)
print(g1)
dev.off()

