setwd("~/git_repos/pseudogene_quantification/data/grasshopper/")
sampDat <- read.csv("sampleData.csv", sep='\t')
str(sampDat)
plot(sampDat$Long, sampDat$Lat, col = sampDat$Karyotype, pch=as.numeric(sampDat$Mitotype))

# available from https://drive.google.com/drive/folders/17dnXkQKlF_fcqqETrHco5cVfF3R7kty0
alt <- read.table("~/Dropbox/elevations/srtm_38_04.asc.gz", sep = " ", skip = 6)
alt <- as.matrix(alt)

library(rgl)

laMin <- 44.2
laMax <- 44.6
loMin <- 6.1
loMax <- 6.8

xllcorner <- 4.9995835028038
yllcorner <- 39.999583575447
cellsize <- 0.00083333333333333


xmin <- floor((loMin - xllcorner)/ cellsize)
xmax <- ceiling((loMax - xllcorner)/ cellsize)
ymin <- abs(floor((laMin - yllcorner)/ cellsize) - 6001)
ymax <- abs(ceiling((laMax - yllcorner)/ cellsize) - 6001)

area <- t(alt)[xmin:xmax, ymin:ymax]
image(area, col=terrain.colors(20))
contour(area, add=T, levels = c(1000,1500,2000))

laConvert <- function(x){
  (x-laMin)/(laMax-laMin)
}
loConvert <- function(x){
  (x-loMin)/(loMax-loMin)
}

laConvert(seq(laMin, laMax, 0.1))
pdf("PodiMap.pdf")
image(area, col=grey.colors(20)[3:20],
      xaxt = 'n',
      yaxt = 'n', 
      xlab="Longitude",
      ylab="Latitude",
      asp=(laMax-laMin)/(loMax-loMin),
      frame = F)

contour(area, add=T, levels = c(1000,1500,2000), col = "#333333")
axis(1, at = loConvert(seq(loMin, loMax, 0.1)), labels = seq(loMin, loMax, 0.1)) # at goes from 0 to 1 in image
axis(2, at = laConvert(seq(laMin, laMax, 0.1)), labels = seq(laMin, laMax, 0.1)) # at goes from 0 to 1 in image
combs <- interaction(sampDat$Mitotype, sampDat$Karyotype)
points(loConvert(sampDat$Long), laConvert(sampDat$Lat), pch = as.numeric(combs), lwd=2)
legend("topright", pch=1:6, legend = levels(combs))
dev.off()
table(sampDat$Karyotype)

