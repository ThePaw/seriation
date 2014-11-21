#!/usr/bin/R
# Bar Plots

setwd("./results/csv")

## R100
# Hours

y=read.csv("R100Hours.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f

pdf("R100Hours.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Duration of 100 iterations in hours, n=100", xlab="", xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rho

y=read.csv("R100Rho.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R100Rho.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Average Spearman Rho, n=100", xlab="", ylim=c(0, 1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rho2

y=read.csv("R100Rho2.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R100Rho2.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Average Spearman Rho, n=100", xlab="", ylim=c(0.90, 1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rank Hits

y=read.csv("R100RankHits.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R100RankHits.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Proportion of Rank Hits, n=100", xlab="", ylim=c(0,1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()



# Hits

y=read.csv("R100Hits.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R100Hits.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Proportion of Hits, n=100", xlab="", ylim=c(0,0.25), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

setwd("../..")

