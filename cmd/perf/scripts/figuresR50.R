#!/usr/bin/R
# Bar Plots

setwd("./results/csv")

## 50
# Sec

y=read.csv("R50Sec.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f

pdf("R50Sec.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Duration of 100 iterations in seconds, n=50", xlab="", xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Sec2

y=read.csv("R50Sec2.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f

pdf("R50Sec2.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Duration of 100 iterations in seconds, n=50", xlab="", xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rho

y=read.csv("R50Rho.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R50Rho.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Average Spearman Rho, n=50", xlab="", ylim=c(0, 1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rho2

y=read.csv("R50Rho2.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R50Rho2.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Average Spearman Rho, n=50", xlab="", ylim=c(0.95, 1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

# Rank Hits

y=read.csv("R50RankHits.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R50RankHits.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Proportion of Rank Hits, n=50", xlab="", ylim=c(0,1), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()



# Hits

y=read.csv("R50Hits.csv", header=FALSE)
f=t(as.vector(y$V1))
x=t(as.matrix(y$V2))
colnames(x) <- f
x

pdf("R50Hits.pdf")
par(pin=c(5, 2.5)) 
par(las=2) 
barplot(x,   main="Proportion of Hits, n=50", xlab="", ylim=c(0,0.8), xpd=FALSE, cex.names=0.8, horiz=FALSE, space=0.5)  
dev.off()

setwd("../..")

