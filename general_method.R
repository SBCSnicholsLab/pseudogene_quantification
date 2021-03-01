################################################################################
# A method to estimate the nuclear frequency of exogenous pseudogenes from     #
#   population-level low-pass sequencing data                                  #
################################################################################

library(lme4) 

# Read in data ####
mainDF <- read.table("transformedData.csv")
head(mainDF)


# The SNPs will have to be sub-sampled for analysis. Generate weights,
#  biasing towards high freq and discrading very low-freq ones.

wweights <- tapply(mainDF$AltProp, mainDF$Position, mean)
hist(wweights, breaks=300)
abline(v=0.01)
wweights[wweights < 0.01] <- 0

# Weights and 

set.seed(12345)
nloci = 200
loci <- sample(names(wweights), nloci, prob = wweights)

subDF <- mainDF[mainDF$Position %in% loci, ]
head(subDF)

# pre-analysis to detect aberrant SNPs
mod4 <- lmer(ylog ~ Position + xnqlogis + xnqlogis:Position + (1 | Sample) ,
             offset = xnqlogis,
             data = subDF)

# Try linear model if model matrix is rank deficient
mod4lm <- lm(ylog ~ Position + xnqlogis + xnqlogis:Position ,
             offset = xnqlogis,
             data = subDF)

#plot(mod4)
longSummary <- summary(mod4)
longSummary <- summary(mod4lm)
interactionTerms <- grep(":xnqlogis", rownames(longSummary$coefficients))
slopeFitted <- longSummary$coefficients[which("xnqlogis" == rownames(longSummary$coefficients)), 1]
hist(longSummary$coefficients[interactionTerms, 1] + slopeFitted, breaks = 25)
slopes <- longSummary$coefficients[interactionTerms, 1] + slopeFitted

# rogues #####

rogueColours <- rep("#00000000", nloci)
rogueColours[which(slopes > 2)] <- 2
rogueColours[which(slopes < 2)] <- 1
#rogueSNPs <- which(abs(slopes) > 0.5)
rogueSNPs <- which(abs(slopes) > 1.5)
rogueSNPs <- as.character(subDF$Position)[rogueSNPs]
goodDat <- subDF[!(subDF$Position %in% rogueSNPs),]
goodDat$Position <- as.character(goodDat$Position)
unique(goodDat$Sample)
goodDat <- goodDat[!(is.na(goodDat$ylog) | is.na(goodDat$xnqlogis)), ]
goodDatTable <- table(goodDat$Position)

hist(goodDatTable)
ltt <- names(goodDatTable[goodDatTable < 2])
goodDat <- goodDat[!goodDat$Position %in% ltt, ]
# goodDat <- subDF


lmod5 <- lmer(ylog ~ 0 + Position + (1 | Sample),
              offset = xnqlogis,
              data = goodDat)


intercepts5 <- summary(lmod5)$coefficients[,1]
summary(intercepts5)
#hist(intercepts5)
iranks <- rank(intercepts5, ties.method = "random")
dim(goodDat)


with(goodDat,{
  plot(ylog~xnqlogis, col=rep(topo.colors(nloci)[iranks], nSamp), pch=1, xlim=c(0,10), ylim=c(-10,0),
       xlab="Odds ratio of data mapped",
       ylab="Allele frequency",
       main="",
       xaxt = 'n',
       yaxt = 'n'
  )
  
  axis(1, at = log(2^(0:7*2)), labels = 2^(0:7*2)) 
  abline(v = log(2^(0:7*2)), col = "grey", lty=3)
  ll <- expression("0", "10"^-1, "10"^-3,"10"^-3,"10"^-4,"10"^-5,"10"^-6)
  axis(2, at = log(10^(0:-6)), labels = ll)
  abline(h = log(10^(0:-6)), col = "grey", lty=3)
  abline(h=max(intercepts5))
  abline(v=max(xnqlogis))
  abline(max(intercepts5), 1,  lty=2)
  text(c(2, 7), c(0,0), c(
       paste0("Lower-bound: ", round(exp(max(intercepts5)), digits = 6)),
       paste0("Upper-bound: ", round(plogis(-max(goodDat$xnqlogis)), digits = 6))
       )
       
  )
  
  sapply(1:length(intercepts5), function(x){
    abline(intercepts5[x], 1, col = topo.colors(nloci, alpha = 0.2)[iranks][x])
  })
})




# Richard's suggestion
with(goodNewDat,{
  plot(ylog~xnqlogis, col=rep(topo.colors(nloci)[iranks], 46), pch=1,
       #xlim=c(0,8), ylim=c(-8,0),
       xlab="(Unmapped reads / Mapped)",
       ylab="Allele Frequency (log scale)",
       main="Estimating nuMt Proportions",
       cex = 0.2,
       xlim = c(0, max(xnqlogis) * 1.01),
       ylim = c(min(intercepts5), max(ylog, na.rm = T) * 1.1),
       xaxt = 'n',
       yaxt = 'n'
  )
  axis(1, at = log(2^(0:3*3)), labels = 2^(0:3*3)) 
  
  ll <- expression("0", "10"^-1, "10"^-3,"10"^-3,"10"^-4,"10"^-5,"10"^-6)
  axis(2, at = log(10^(0:-6)), labels = ll)
  
  sapply(1:length(intercepts5), function(x){
    abline(intercepts5[x], 1, col = topo.colors(nloci, alpha = 0.2)[iranks][x])
  })
  abline(v = 0, lty = 3)
  #plot(ylog~xlog, col=rep(topo.colors(nloci)[iranks], 46), pch=1)
})




