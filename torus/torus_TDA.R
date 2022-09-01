library("transport")
library(igraph)
library(devtools)
library(dimension)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## ================== persistence diagrams

tda_Z_H0 = read.table("tda_Z_H0.dat")
tda_Z_H1 = read.table("tda_Z_H1.dat")

tda_PC_d3_H0 = read.table("tda_PC_d3_H0.dat")
tda_PC_d3_H1 = read.table("tda_PC_d3_H1.dat")

tda_PC_d20_H0 = read.table("tda_PC_d20_H0.dat")
tda_PC_d20_H1 = read.table("tda_PC_d20_H1.dat")

tda_Z_H0['col'] = 'black'
tda_Z_H1['col'] = 'red'

tda_PC_d3_H0['col'] = 'black'
tda_PC_d3_H1['col'] = 'red'

tda_PC_d20_H0['col'] = 'black'
tda_PC_d20_H1['col'] = 'red'

tda_Z = rbind(tda_Z_H0,tda_Z_H1)
tda_PC_d3 = rbind(tda_PC_d3_H0,tda_PC_d3_H1)
tda_PC_d20 = rbind(tda_PC_d20_H0,tda_PC_d20_H1)

tda_Z$V2 <- ifelse(is.infinite(tda_Z$V2), 1.5, tda_Z$V2)
tda_PC_d3$V2 <- ifelse(is.infinite(tda_PC_d3$V2), 1.5, tda_PC_d3$V2)
tda_PC_d20$V2 <- ifelse(is.infinite(tda_PC_d20$V2), 1.5, tda_PC_d20$V2)

dev.new(width=10, height=3.5)
# pdf(file = 'torus_tda.pdf',width = 10,height = 3.5)
par(mfrow=c(1,3))
plot(tda_Z[,1:2], col=tda_Z$col, pch=16, xlim=range(c(tda_Z[,1:2])), ylim=range(c(tda_Z[,1:2])), cex=.5,
     xlab = 'Birth', ylab = 'Death')
abline(a = 0, b=1, lty=1)
legend("bottomright", legend=c("Dim 0", "Dim 1"), pch=c(16,16), col=c("black", "red"))
# abline(h=  1.5, lty=1)

plot(tda_PC_d3[,1:2], col=tda_PC_d3$col, pch=16, xlim=range(c(tda_PC_d3[,1:2])), ylim=range(c(tda_PC_d3[,1:2])), cex=.5,
     xlab = 'Birth', ylab = 'Death')
abline(a = 0, b=1, lty=1)
legend("bottomright", legend=c("Dim 0", "Dim 1"), pch=c(16,16), col=c("black", "red"))

plot(tda_PC_d20[,1:2], col=tda_PC_d20$col, pch=16, xlim=range(c(tda_PC_d20[,1:2])), ylim=range(c(tda_PC_d20[,1:2])), cex=.5,
     xlab = 'Birth', ylab = 'Death')
abline(a = 0, b=1, lty=1)
legend("bottomright", legend=c("Dim 0", "Dim 1"), pch=c(16,16), col=c("black", "red"))

# dev.off()
