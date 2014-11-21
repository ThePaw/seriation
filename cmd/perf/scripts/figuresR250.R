#!/usr/bin/R
# Bar Plots

setwd("./results/csv")

## R250
# Hours

y=read.csv("R250Hours.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f

pdf("R250Hours.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Duration of 100 iterations in hours, n=250", xlab="", xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rho

y=read.csv("R250Rho.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R250Rho.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Average Spearman Rho, n=250", xlab="", ylim=c(0, 1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rho2

y=read.csv("R250Rho2.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R250Rho2.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Average Spearman Rho, n=250", xlab="", ylim=c(0.90, 1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()


setwd("../..")

