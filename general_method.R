################################################################################
# A method to estimate the nuclear frequency of exogenous pseudogenes from     #
#   population-level low-pass sequencing data                                  #
################################################################################

library(lme4) 


# The analysis function
extrasIns <- function(x, seed, main=""){
  # The SNPs will have to be sub-sampled for analysis. Generate weights,
  #  biasing towards high freq and discrading very low-freq SNPs.
  wweights <- tapply(x$AltProp, x$Position, mean)
  # hist(wweights, breaks=300, xlim=c(0,0.1))
  # abline(v=0.01)
  wweights[wweights < 0.01] <- 0
  
  # Weights
  if(hasArg(seed)) set.seed(seed)
  nloci = 200
  loci <- sample(names(wweights), nloci, prob = wweights)
  
  subDF <- x[x$Position %in% loci, ]
  
  
  # pre-analysis to detect aberrant SNPs
  # mod4 <- lmer(ylog ~ Position + xnqlogis + xnqlogis:Position + (1 | Sample) ,
  #              offset = xnqlogis,
  #              data = subDF)
  
  # Try linear model if model matrix is rank deficient
  mod4lm <- lm(ylog ~ Position + xnqlogis + xnqlogis:Position ,
               offset = xnqlogis,
               data = subDF)
  
  #plot(mod4)
  # longSummary <- summary(mod4)
  longSummary <- summary(mod4lm)
  interactionTerms <- grep(":xnqlogis", rownames(longSummary$coefficients))
  slopeFitted <- longSummary$coefficients[which("xnqlogis" == rownames(longSummary$coefficients)), 1]
  slopes <- longSummary$coefficients[interactionTerms, 1] + slopeFitted
  
  # rogues #####
  
  # hist(slopes, breaks = 100)
  # hist(slopes, breaks = 1000, xlim=c(-10,10))
  bpsSlopes <- boxplot.stats(slopes)$stats
  #abline(v=bpsSlopes)
  rogueSNPs <- which(slopes > bpsSlopes[5] | slopes < bpsSlopes[1])
  rogueSNPs <- as.character(subDF$Position)[rogueSNPs]
  goodDat <- subDF[!(subDF$Position %in% rogueSNPs),]
  goodDat$Position <- as.character(goodDat$Position)
  
  goodDat <- goodDat[!(is.na(goodDat$ylog) | is.na(goodDat$xnqlogis)), ]
  goodDat <- goodDat[!goodDat$AltProp > 0.5, ]
  
  goodDatTable <- table(goodDat$Position)
  # hist(goodDatTable)
  
  ltt <- names(goodDatTable[goodDatTable < 5])
  goodDat <- goodDat[!goodDat$Position %in% ltt, ]
  
  # goodDat <- subDF
  
  
  lmod5 <- lmer(ylog ~ 0 + Position + (1 | Sample),
                offset = xnqlogis,
                data = goodDat)
  
  
  intercepts5 <- summary(lmod5)$coefficients[,1]
  summary(intercepts5)
  #hist(intercepts5)
  
  
  # rank intercepts, so allele freqs and slopes can be coloured in
  iranks <- rank(intercepts5, ties.method = "random")
  rankDF <- data.frame(Position=substr(names(iranks), 9, 20), rank=iranks)
  goodDat <- merge(goodDat, rankDF, by="Position")
  
  
  
  with(goodDat,{
    nSamp = length(unique(goodDat$Sample))
    remainingLoci = length(unique(goodDat$Position))
    plot(ylog~xnqlogis, col=rainbow(remainingLoci)[rank], pch=1, xlim=c(0,10), ylim=c(-10,0),
         xlab="(Un-mapped data / mapped)",
         ylab="Allele frequency",
         main=main,
         xaxt = 'n',
         yaxt = 'n'
    )

    axis(1, at = log(2^(0:7*2)), labels = 2^(0:7*2))
    abline(v = log(2^(0:7*2)), col = "grey", lty=3)
    ll <- expression("0", "10"^-1, "10"^-2, "10"^-3, "10"^-4, "10"^-5, "10"^-6)
    axis(2, at = log(10^(0:-6)), labels = ll)
    abline(h = log(10^(0:-6)), col = "grey", lty=3)
    abline(h=max(intercepts5))
    abline(v=max(xnqlogis))
    abline(max(intercepts5), 1,  lty=2)
    estLB <<- exp(max(intercepts5))/ (1+exp(max(intercepts5)))
    estUB <<- plogis(-max(goodDat$xnqlogis))
    text(c(2, 7), c(0,0), c(
      paste0("Lower-bound: ", round(estLB*100, digits = 3), "%"),
      paste0("Upper-bound: ", round(estUB*100, digits = 3), "%")
    )

    )

    sapply(1:length(intercepts5), function(x){
      abline(intercepts5[x], 1, col = rainbow(remainingLoci, alpha = 0.2)[iranks][x])
    })
    print(paste0("Estimate based on ", length(intercepts5), " loci after filtering."))
    #print(paste0("Lower-bound: ", estLB))
    #print(paste0("Upper-bound: ", estUB))
  })
  return(c(lower=estLB, upper=estUB))
}


# Get data online ####
# human
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/human/transformedData.csv",
              destfile = "human.csv")
humanDF <- read.table("human.csv", stringsAsFactors = F)

#hopper
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/grasshopper/transformedData.csv",
              destfile = "hopper.csv")
hopperDF <- read.table("hopper.csv", stringsAsFactors = F)

#parrot
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/parrot/transformedData.csv",
              destfile = "parrot.csv")
parrotDF <- read.table("parrot.csv", stringsAsFactors = F)

# Alternatively, read in data from repo ####
# setwd("~/git_repos/pseudogene_quantification/")
# # human data set
# humanDF <- read.table("data/human/transformedData.csv", stringsAsFactors = F)
# # parrot data set
# parrotDF <- read.table("data/parrot/transformedData.csv", stringsAsFactors = F)
# # grasshopper data set
# hopperDF <- read.table("data/grasshopper/transformedData.csv", stringsAsFactors = F)


extrasIns(humanDF, seed = 12345, main="Human")

extrasIns(hopperDF, seed = 12345, main="Grasshopper")

extrasIns(parrotDF, seed = 12345, main="Parrot")
