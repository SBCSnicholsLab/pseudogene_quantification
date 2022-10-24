
setwd("~/Dropbox/Richard_share/mtAnalysis2020/@manuscript/Figures/")
sampDat <- read.csv("~/git_repos/pseudogene_quantification/data/grasshopper/sampleData.csv", sep='\t')
str(sampDat)
plot(sampDat$Long, sampDat$Lat, pch=as.numeric(as.factor(sampDat$Mitotype)))

# available from https://drive.google.com/drive/folders/17dnXkQKlF_fcqqETrHco5cVfF3R7kty0
alt <- read.table("~/Dropbox/elevations/srtm_38_04.asc.gz", sep = " ", skip = 6)
alt <- as.matrix(alt)

#library(rgl)

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

image(area, col=c("#777777", "#888888", "#AAAAAA", "#CCCCCC"), breaks=c(0,1000,1500, 2000, 5000))
contour(area, add=T, levels = c(1000,1500,2000))

laConvert <- function(x){
  (x-laMin)/(laMax-laMin)
}
loConvert <- function(x){
  (x-loMin)/(loMax-loMin)
}

laConvert(seq(laMin, laMax, 0.1))
#pdf("PodiMap.pdf")
image(area, col=grey.colors(20)[3:20],
      xaxt = 'n',
      yaxt = 'n',
      xlab="Longitude",
      ylab="Latitude",
      asp=(laMax-laMin)/(loMax-loMin),
      frame = F)

#png("PodiMapBlank.png", width=4000, height=4000)
image(area, col=c("#777777", "#888888", "#AAAAAA", "#CCCCCC"), breaks=c(0,1000,1500, 2000, 5000),
      xaxt = 'n',
      yaxt = 'n', 
      xlab="Longitude",
      ylab="Latitude",
      asp=(laMax-laMin)/(loMax-loMin),
      frame = F)
contour(area, add=T, levels = c(1000,1500,2000), col = "#111111", lwd=10, labcex = 5)
#dev.off()
axis(1, at = loConvert(seq(loMin, loMax, 0.1)), labels = seq(loMin, loMax, 0.1)) # at goes from 0 to 1 in image
axis(2, at = laConvert(seq(laMin, laMax, 0.1)), labels = seq(laMin, laMax, 0.1)) # at goes from 0 to 1 in image

points(loConvert(sampDat$Long), laConvert(sampDat$Lat), pch = as.numeric(combs), lwd=2)
legend("topright", pch=1:6, legend = levels(as.factor(sampDat$Mitotype)))
dev.off()


#pdf("contours.pdf")
plot(c(0,0), c(1,1), type='n', xlim=c(0,1), ylim=c(0,1),asp=(laMax-laMin)/(loMax-loMin))
contour(area, add=T, levels = c(1000,1500,2000), col = "#111111")
#dev.off()
