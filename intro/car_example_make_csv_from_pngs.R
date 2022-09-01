library(png)
library(plotly)
library(RSpectra)

nangles = 72
npixels = 384*288 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data_mat <- matrix(nrow = npixels, ncol = nangles)

##160 is the label of the car images

filenames = sprintf("160_r%d.png", (0:72*5))
for(i in 1:nangles){
 path <- paste0("./car images/", filenames[i])
 data_mat[, i] <- c(readPNG(path))
}
write.csv(data_mat,"car_pixels.csv", row.names = FALSE, col.names=FALSE)


