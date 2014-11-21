#!/usr/bin/R
# Bar Plots

setwd("./results/csv")

## R500
# Hours

y=read.csv("R500Hours.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f

pdf("R500Hours.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Duration of 100 iterations in hours, n=500", xlab="", xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rho

y=read.csv("R500Rho.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R500Rho.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Average Spearman Rho, n=500", xlab="", ylim=c(0, 1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()


setwd("../..")

