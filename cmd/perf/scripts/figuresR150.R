#!/usr/bin/R
# Bar Plots

setwd("./results/csv")

## R150
# Hours

y=read.csv("R150Hours.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f

pdf("R150Hours.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Duration of 100 iterations in hours, n=150", xlab="", xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Hours2

y=read.csv("R150Hours2.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f

pdf("R150Hours2.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Duration of 100 iterations in hours, n=150", xlab="", xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Hours3

y=read.csv("R150Hours3.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f

pdf("R150Hours3.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Duration of 100 iterations in hours, n=150", xlab="", xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()




# Rho

y=read.csv("R150Rho.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R150Rho.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Average Spearman Rho, n=150", xlab="", ylim=c(0, 1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rank Hits

y=read.csv("R150RankHits.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R150RankHits.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Proportion of Rank Hits, n=150", xlab="", ylim=c(0,1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

setwd("../..")

