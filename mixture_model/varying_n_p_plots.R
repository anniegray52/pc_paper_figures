library(scatterplot3d)
library(tidyverse)

source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

## need to run mixture_model.ipynb first 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
ns_df = read.csv(file = 'ns_df_p200.csv')
ps_df = read.csv(file = 'ps_df_n200.csv')

ns = c(100, 200, 500, 1000)
ps = c(100, 500, 1000, 10000)
a = 45

# pdf(file = paste(c(home, 'varying_n.pdf'), collapse="/"),width = 20, height = 5)
dev.new(width=20, height=5)
par(mfrow=c(1,4))
for (i in ns){
  df = ns_df %>% filter(ns_df$n == i)
  s3d <- scatterplot3d(df[, 1:3], pch = "",  angle = a,grid=FALSE, box=FALSE, 
                       xlab= '',ylab= '',zlab= '', main = paste(c('n =', i),collapse=" "), 
                       xlim=c(min(ns_df[,1]),max(ns_df[,1])),
                       ylim=c(min(ns_df[,2]),max(ns_df[,2])),
                       zlim=c(min(ns_df[,3]),max(ns_df[,3])),cex.main=4,tick.marks = FALSE)
  addgrids3d(df[, 1:3], grid = c("xy", "xz", "yz"),
             xlim=c(min(ns_df[,1]),max(ns_df[,1])),
             ylim=c(min(ns_df[,2]),max(ns_df[,2])),
             zlim=c(min(ns_df[,2]),max(ns_df[,3])), angle = a)
  s3d$points3d(df[, 1:3],cex=3, col ='black', bg = df$type,lwd=.1, pch = 21)
}
# dev.off()



# pdf(file = paste(c(home, 'varying_p.pdf'), collapse="/"),width = 20,height = 5)
dev.new(width=20, height=5)
par(mfrow=c(1,4))
for (i in ps){
  df = ps_df %>% filter(ps_df$p == i)
  s3d <- scatterplot3d(df[, 1:3], pch = "",  angle = a,grid=FALSE, box=FALSE, 
                       xlab= '',ylab= '',zlab= '', main = paste(c('p =', i), collapse=" "),
                       xlim=c(min(ps_df[,1]),max(ps_df[,1])),
                       ylim=c(min(ps_df[,2]),1.4),#max(ps_df[,2])),
                       zlim=c(min(ps_df[,3]),max(ps_df[,3])),cex.main=4,tick.marks = FALSE)
  addgrids3d(df[, 1:3], grid = c("xy", "xz", "yz"), 
             xlim=c(min(ps_df[,1]),max(ps_df[,1])),
             ylim=c(min(ps_df[,2]),1.4),#max(ps_df[,2])),
             zlim=c(min(ps_df[,3]),max(ps_df[,3])),angle = a)
  s3d$points3d(df[, 1:3],cex=3 , col ='black', bg = df$type,lwd=.1, pch = 21)
}
# dev.off()

