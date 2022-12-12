


rPar <- rainbowPlot(parrotDF, seed=12345)
# remove multiple outliers
toRemove <- interceptPositionPlot(rPar, highlightOutliers=T)
parrotDF2 <- parrotDF[!parrotDF$Position %in% toRemove,]
rainbowPlot(parrotDF2, seed=12345)


rHum <- rainbowPlot(humanDF, seed=12345)
# automatic outlier detection would remove SNPs that look OK
remHum <- interceptPositionPlot(rHum, highlightOutliers = T)

# select single outlier "by hand"
which.max(coef(summary(rHum$lmer.model))[,1])

humanDF2 <- humanDF[humanDF$Position != 310,] # remove "by hand"
rainbowPlot(humanDF2, seed=12345)
divEst(parrotFX)

# generate plots for manuscript
#pdf("hopper.pdf", width=6, height=7)

#download.file("https://tinyurl.com/4mtrbkzc", destfile = "hopper.csv")
hopperDF <- read.table("hopper.csv")

hopperFit <- rainbowPlot(hopperDF,
                         seed = 12345, printout = FALSE, title = "Grasshopper")
#dev.off()


#pdf("parrot.pdf", width=6, height=7)
rainbowPlot(parrotDF2,
            seed = 12345, printout = FALSE, title = "Parrot")
#dev.off()

#pdf("human.pdf", width=6, height=7)
rainbowPlot(humanDF2,
            seed = 12345, printout = FALSE, title = "Human")
#dev.off()
