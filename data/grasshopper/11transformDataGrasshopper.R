######################################
# Prepare data for general method ####
######################################

################
# mapping depths

# adjust directory
setwd("~/git_repos/pseudogene_quantification/data/grasshopper/")
# setwd("~/Documents/GitHub/pseudogene_quantification/data/grasshopper")

dir()
nSamp = 52
# a utility function to rename the mapping depth vals
getNums <- function(x){
  a <- x[,3]
  names(a) <- sapply(x[,1], function(y){
    strsplit(y, "j")[[1]][1]
  })
  return(a)
}
nMapped <- read.table("bpMapped", sep = "\t", stringsAsFactors = F)
head(nMapped)

nTot <- read.table("bpTotal", sep = "\t", stringsAsFactors = F)
head(nTot)

nMapped <- getNums(nMapped)
nTot <- getNums(nTot)

mappingProp <- nMapped / nTot
hist(mappingProp, breaks=20)


#################
# Gentotype calls

gts <- read.table("genotypes.csv",header = T, check.names = F)
head(gts)
# counts of each allele, some are missing (-1)
table(as.vector(unlist(gts[,1:nSamp])))

# remove lowest coverage samples
gtsHC <- gts # none removed

# visualise missingness
#image(as.matrix(gtsHC[,1:nSamp]) == -1)

# remove any site with missing data
posMissing <- gtsHC$POS[apply(gtsHC, 1, function(x) -1 %in% x)]
gtsHC <- gtsHC[!gtsHC$POS %in% posMissing,]

###############
# Allele counts

allCounts <- read.table("alleleCounts.csv")
head(allCounts)
allCountsHC <- allCounts # opportunity to remove individuals if desired

# remove sites with missing alleles identified above
allCountsHC <- allCountsHC[!allCountsHC$position %in% posMissing,]

# split allele counts into separate DFs, one per allele
gtList <- lapply(1:4, function(x) {
  allCountsHC[allCountsHC$allele==x,]
})

lapply(gtList, dim)
# 7477*4 = 29908

# flatten into four vectors
rm(c) # In case somebody named a variable c, which would shadow the c function.
gtVectors <- lapply(gtList, function(x){
  do.call(c,x[,1:nSamp])
})
#index of samples
sampleIndex <- rep(colnames(gtsHC)[1:nSamp], each=nrow(gtList[[1]]))
str(gtVectors)

# long(er) data frame for allele counts with sample index
allDF <- data.frame(sample=sampleIndex, pos = paste0("S",rep(sprintf("%05d",gtsHC$POS), nSamp)), 
                    ref=gtVectors[[1]], alt=gtVectors[[2]] + gtVectors[[3]] + gtVectors[[4]])
head(allDF)
str(allDF)

# Depth of the alt allels (summed)
siteDeps <- rowSums(allDF[,3:4])
altRatio <- allDF$alt / siteDeps



# DF that wih all releavant data for the anaysis:
mainDF <- data.frame(Sample = as.character(allDF$sample),
                     Position=allDF$pos,
                     AltProp=altRatio,
                     stringsAsFactors = F, row.names = NULL)

# Add mapping rate, make a DF so we can use "merge"
mappingProp
propDF <- data.frame(Sample = names(mappingProp), mappingrate=mappingProp, stringsAsFactors = F)
head(propDF)
mainDF <- merge(mainDF, propDF, by="Sample")

head(mainDF)
mainDF$ylog <- log(mainDF$AltProp)
mainDF$ylog[is.infinite(mainDF$ylog)] <- NA # replace -Inf by NA, makes plotting easier
mainDF$xnqlogis <- -qlogis(mainDF$mappingrate)
mainDF$xlog <- log(mainDF$mappingrate)
head(mainDF)

# remove NA lines

mainDF <- mainDF[!is.na(mainDF$AltProp),]
mainDF <- mainDF[!is.na(mainDF$ylog),]
mainDF <- mainDF[!is.na(mainDF$xnqlogis),]

# Write out data
#write.table(mainDF, "transformedData.csv")



#################################################
# Diverged populations and fixed differences ####
#################################################


# population divergence ####


dim(gtsHC)
head(gtsHC)
pc01 <- prcomp(t(gtsHC[,1:nSamp]))

# PC1 accounts for >60% of the variance, PC2 for <3%
summary(pc01)
str(pc01)
pc01$x[,1:2]

# PCA separates two clusters of individuals
plot(pc01$x[,1:2], asp = 6/64, main="PCA mitotypes")

# # optional
# library(ggplot2)
# library(ggrepel) # too many overlapping points
# plotDF <- data.frame(pc01$x[,1:2])
# ggplot(plotDF, aes(PC1, PC2, label=rownames(plotDF))) +
#   geom_text_repel() +
#   geom_point()

# Which sites have fixed differences between the groups? ####

# Get sample names of PCA rgoups
pop0 <- rownames(pc01$x)[pc01$x[,1]>0]
pop1 <- rownames(pc01$x)[pc01$x[,1]<0]

# Get indices
pop0ind <- which(names(gtsHC) %in% pop0)
pop1ind <- which(names(gtsHC) %in% pop1)

# A function to check whether there is sharing of alleles
#  (for an individual locus)
is.fixed.diff <- function(x, ind0, ind1){
  !any(x[ind0] %in% x[ind1])
}

# The first locus in gtsHC is not fixed:
is.fixed.diff(gtsHC[1,], pop1ind, pop0ind)

cov15 <- c(names(which(mappingProp > 0.002)), "POS")


names(gtsHC) %in% cov15
gtsCov15 <- gtsHC[, names(gtsHC) %in% cov15]

pop0.15 <- pop0[pop0 %in% names(gtsCov15)]
pop1.15 <- pop1[pop1 %in% names(gtsCov15)]
pop0ind.15 <- which(names(gtsCov15) %in% pop0.15)
pop1ind.15 <- which(names(gtsCov15) %in% pop1.15)


head(gtsCov15)
fixedInd <- apply(gtsCov15, 1, function(x) is.fixed.diff(x, pop0ind.15, pop1ind.15))
sum(fixedInd)

fixedPos <- paste0("S", sprintf("%05d", gtsHC$POS[fixedInd]))

###################################################
# Additional estimate for diverged populations ####
###################################################

allCountsLonger <- data.frame(sample=sampleIndex, pos = paste0("S",rep(sprintf("%05d",gtsHC$POS), nSamp)), 
           g1=gtVectors[[1]], g2=gtVectors[[2]], g3=gtVectors[[3]], g4=gtVectors[[4]])


head(allCountsLonger)
str(allCountsLonger)

# A matrix to indicate which allele is has the highest count in each individual at each locus (position).
# There are ties, which we will resolve later, comparing across individuals within each sub-population.
maxMat <- t(apply(as.matrix(allCountsLonger[, 3:6]), 1, function(x) (x == max(x))))
dim(maxMat)

allCountsLonger <- data.frame(allCountsLonger, maxMat)
head(allCountsLonger)

# use a join to get the subset of allDF with fixedDifferences
allCountsFixed <- merge(allCountsLonger, data.frame(pos = fixedPos), by="pos")
dim(allCountsLonger)
dim(allCountsFixed)
head(allCountsFixed)

# Adding a column for N (the non-mitolike reads)
# same order?
all(names(nMapped) == names(nTot))
# yes.
# another join to get the numbers of non-mitolike reads
allCountsFixedN <- merge(allCountsFixed, data.frame(sample=names(nMapped),
                                                    N=unname(nTot-nMapped),
                                                    M=unname(nMapped)), by="sample")
dim(allCountsFixed)
dim(allCountsFixedN) # we have not lost any lines

allCountsFixedN$pop <- ifelse(allCountsFixedN$sample %in% pop1, "A", "B")
head(allCountsFixedN)
str(allCountsFixedN)
# get rid of unneeded factor levels
allCountsFixedN$pop <- as.character(allCountsFixedN$pop)
allCountsFixedN$pos <- as.character(allCountsFixedN$pos)
majAlleles <- tapply(1:nrow(allCountsFixedN), list(allCountsFixedN$pos, allCountsFixedN$pop), function(x) {
  a <- colSums(allCountsFixedN[x,7:10])
  which(a == max(a))
})
majAlleles <- data.frame(majAlleles, pos = rownames(majAlleles))
# All major alleles are different between the sub-pops (no 0 in the histogram):
hist(majAlleles$A - majAlleles$B)

allCountsFixedN <- merge(allCountsFixedN, majAlleles, by = "pos")
allCountsFixedN <- allCountsFixedN[,-c(7:10)]
# write.table(allCountsFixedN, "hopperFixed.csv")
# estimation is now done with the divEst function!

# allCountsFixedN <- read.table("hopperFixed.csv")
# 
# head(allCountsFixedN)
# 
# allCountsFixedN$A1 <- allCountsFixedN$A == 1
# allCountsFixedN$A2 <- allCountsFixedN$A == 2
# allCountsFixedN$A3 <- allCountsFixedN$A == 3
# allCountsFixedN$A4 <- allCountsFixedN$A == 4
# allCountsFixedN$B1 <- allCountsFixedN$B == 1
# allCountsFixedN$B2 <- allCountsFixedN$B == 2
# allCountsFixedN$B3 <- allCountsFixedN$B == 3
# allCountsFixedN$B4 <- allCountsFixedN$B == 4
# allCountsFixedN$allDep <- rowSums(allCountsFixedN[,3:6])
# allCountsFixedN$aMaj <- rowSums(allCountsFixedN[,3:6] * allCountsFixedN[,12:15])
# allCountsFixedN$aAlt <- rowSums(allCountsFixedN[,3:6] * !allCountsFixedN[,12:15])
# allCountsFixedN$bMaj <- rowSums(allCountsFixedN[,3:6] * allCountsFixedN[,16:19])
# allCountsFixedN$bAlt <- rowSums(allCountsFixedN[,3:6] * !allCountsFixedN[,16:19])
# # Richard's scores ###################################################
# #####################################################################
# 
# # extract site names & number
# 
# siteNames <- unique(allCountsFixedN$pos)
# nSites <- length(siteNames)
# Ascores <- Bscores <- rep(0,nSites)
# 
# for (i in 1:nSites){# create matrices with counts for the A populations and B populations
#   gmatA <- as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')[,3:6])
#   gmatB <- as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')[,3:6])
#   
#   # create matrices identifying the non mito alleles in each population 
#   # and the mito allele in the opposite population (where it will be non-mito)
#   # For the A population:
#   notMitoAlleleA <- !as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')[,12:15])
#   AmitoAlleleInB <-  as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')[,12:15])
#   # For the B populations
#   notMitoAlleleB <- !as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')[,16:19]) 
#   BmitoAlleleInA <-  as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')[,16:19])
#   mA <- subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')$M
#   mB <- subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')$M
#   nA <- subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')$N
#   nB <- subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')$N
#   
#   Ascores[i] <- sum(rowSums(gmatA * notMitoAlleleA) / rowSums(gmatA) * mA) / sum(nA) +
#                 sum(rowSums(gmatB * AmitoAlleleInB) / rowSums(gmatB) * mB) / sum(nB)  
#   Bscores[i] <- sum(rowSums(gmatB * notMitoAlleleB) / rowSums(gmatB) * mB) / sum(nB) +
#                 sum(rowSums(gmatA * BmitoAlleleInA) / rowSums(gmatA) * mA) / sum(nA)  
#   }
# 
# 
# plot(Ascores, Bscores)
# abline(0,1)
# abline(coef(lm(Bscores ~ Ascores)), col = "Blue")
# abline(v = mean(Ascores), col='Orange')
# abline(h = mean(Bscores), col='Orange')
# 
# noquote(paste("Mean of A Scores: ",
#       signif(mean(Ascores), 3),
#       " (SE ",
#       signif(sqrt(var(Ascores)/nSites),2),
#       ")"
#       ))
# 
# noquote(paste("Mean of B Scores: ",
#               signif(mean(Bscores), 3),
#               " (SE ",
#               signif(sqrt(var(Bscores)/nSites),2),
#               ")"
# ))
# 
# plot(as.numeric(substr(siteNames, 2, 6)), Ascores, type="l",
#      xlab="Position in mitogenome",
#      ylab="Score")
# points(as.numeric(substr(siteNames, 2, 6)), Bscores, col="grey", type="l")
# legend("top", lty=1, col=c("black","grey"), legend=c("Ascores","Bscores"))
# 
# plot(Ascores, Bscores)
# lm0 <- lm(Bscores~Ascores)
# lm1 <- lm(Bscores~Ascores, offset = Ascores)
# summary(lm0)
# summary(lm1)
# abline(lm0)
# abline(0, 1, lty = 2)
# 
# 
# mappingDepths <- data.frame(nMapped, pop=rep("A", length(nMapped)), stringsAsFactors = F)
# mappingDepths$pop[pop0ind]
# mappingDepths$pop[pop1ind] <- "B"
# summary(lm(nMapped ~ pop, data=mappingDepths))
# #plot(lm(nMapped ~ pop, data=mappingDepths))
# 
# hist(Ascores, col = "#0000FF40", main = "Grasshopper", xlab="Scores", xlim=c(0, 1/1000))
# abline(v=mean(Ascores), lwd=2)
# abline(v=mean(Ascores) + c(-1, 1) * sqrt(var(Ascores)/nSites), lty=2)
# # intervals overlap
# # hist(Bscores, add = T, col = "#FF000040")
# # abline(v=mean(Bscores), lwd=2)
# # abline(v=mean(Bscores) + c(-1, 1) * sqrt(var(Ascores)/nSites), lty=3)
# # legend("topright", fill=c("#0000FF40","#FF000040"), legend=c("A","B"))
