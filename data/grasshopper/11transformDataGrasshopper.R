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

################################################
# New statistic ####
################################################

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
# update the major allele matrix
oneTrue <- function(x){
  a = c(F, F, F, F)
  a[x] <- T
  a
}
# oneTrue(2)
# This is only 6k rows, so we can use sapply.
majMat <- as.matrix(t(sapply(1:nrow(allCountsFixedN), function(x) {
  c(oneTrue(allCountsFixedN[x, ]$A),
    oneTrue(allCountsFixedN[x, ]$B)
  )})))
# only one major allele now?
all(rowSums(majMat) == 2)
# yes, there are two because there's one major allele in A and one in B
colnames(majMat) <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4")
dim(majMat)
dim(allCountsFixedN)
head(allCountsFixedN)
allCountsFixedN <- data.frame(allCountsFixedN[, -c(7:10)], majMat)
head(allCountsFixedN)
allCountsFixedN$allDep <- rowSums(allCountsFixedN[,3:6])
head(allCountsFixedN)
# get numbers for each locus and allele aMaj and bAlt are sometimes indentical,
# but only if there were only two alleles.
aMaj <- rowSums(allCountsFixedN[,3:6] * allCountsFixedN[,12:15])
bMaj <- rowSums(allCountsFixedN[,3:6] * allCountsFixedN[,16:19])
aAlt <- rowSums(allCountsFixedN[,3:6] * !allCountsFixedN[,12:15])
bAlt <- rowSums(allCountsFixedN[,3:6] * !allCountsFixedN[,16:19])
allCountsFixedN <- data.frame(allCountsFixedN, aMaj, aAlt, bMaj, bAlt)
head(allCountsFixedN)
allCountsFixedN[allCountsFixedN$pop == "A" & allCountsFixedN$pos == "S00535",-1]
allCountsFixedN[allCountsFixedN$pop == "A" & allCountsFixedN$pos == "S01184",-1]
str(allCountsFixedN)
# testSet <- allCountsFixedN[allCountsFixedN$pos == "S00535", ]
# tapply(1:nrow(testSet), testSet$pop, function(y){
#   colSums(testSet[y,c(7,8, 20:24)])
# })

# newStat <- do.call(rbind, 
#                    tapply(1:nrow(allCountsFixedN), allCountsFixedN$pos, function(x){
#   a <- allCountsFixedN[x,]
#   sums <- tapply(1:nrow(a), a$pop, function(y){
#     colSums(a[y,c(7, 19:22)])
#   })
#   stats <- (c(aAlt = unname(sums$A["aAlt"]/sums$A["N"] + sums$B["bMaj"]/sums$B["N"]),
#               bAlt = unname(sums$B["bAlt"]/sums$B["N"] + sums$A["aMaj"]/sums$A["N"])
#   ))
#   
# }))
newStat <- do.call(rbind, 
                   tapply(1:nrow(allCountsFixedN), allCountsFixedN$pos, function(x){
                     a <- allCountsFixedN[x,]
                     sums <- tapply(1:nrow(a), a$pop, function(y){
                       colSums(a[y,c(7,8, 20:24)])
                     })
                     stats <- (c(aAlt = unname(sums$A["aAlt"]/sums$A["allDep"]*sums$A["M"]/sums$A["N"] + sums$B["bMaj"]/sums$B["allDep"]*sums$B["M"]/sums$B["N"]),
                                 bAlt = unname(sums$B["bAlt"]/sums$B["allDep"]*sums$B["M"]/sums$B["N"] + sums$A["aMaj"]/sums$A["allDep"]*sums$A["M"]/sums$A["N"])
                     ))
                     
                   }))



head(newStat)
# plot 
plot(newStat)
grid()

# log axes
plot(newStat, log="xy")
grid()

# Histograms
hist(newStat)
hist(newStat[,1], col = "#FF000040", add=T)
hist(newStat[,2], col = "#00FF0040", add=T)
legend("topright", fill=c("white", "#FF000040", "#00FF0040"),
       legend = c("all", "aAlt", "bAlt"))


# Richards attempt at the new statistics ###################################################
############################################################################################

# extract site names & number
siteNames <- unique(allCountsFixedN$pos)
nSites <- length(siteNames)
Ascores <- Bscores <- rep(0,nSites)

for (i in 1:nSites){# create matrices with counts for the A populations and B populations
  gmatA <- as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')[,3:6])
  gmatB <- as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')[,3:6])
  
  # create matrixes identifying the non mito alleles in each population 
  # and the mito allele in the opposite population (where it will be non-mito)
  # For the A population:
  notMitoAlleleA <- !as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')[,12:15])
  AmitoAlleleInB <-  as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')[,12:15])
  # For the B populations
  notMitoAlleleB <- !as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')[,16:19]) 
  BmitoAlleleInA <-  as.matrix(subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')[,16:19])
  mA <- subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')$M
  mB <- subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')$M
  nA <- subset(allCountsFixedN, pos==siteNames[i] & pop == 'A')$N
  nB <- subset(allCountsFixedN, pos==siteNames[i] & pop == 'B')$N
  
  Ascores[i] <- sum(rowSums(gmatA * notMitoAlleleA) / rowSums(gmatA) * mA) / sum(nA) +
                sum(rowSums(gmatB * AmitoAlleleInB) / rowSums(gmatB) * mB) / sum(nB)  
  Bscores[i] <- sum(rowSums(gmatB * notMitoAlleleB) / rowSums(gmatB) * mB) / sum(nB) +
                sum(rowSums(gmatA * BmitoAlleleInA) / rowSums(gmatA) * mA) / sum(nA)  
  }


plot(Ascores, Bscores)
abline(0,1)
abline(coef(lm(Bscores ~ Ascores)), col = "Blue")
abline(v = mean(Ascores), col='Orange')
abline(h = mean(Bscores), col='Orange')

noquote(paste("Mean of A Scores: ",
      signif(mean(Ascores), 3),
      " (SE ",
      signif(sqrt(var(Ascores)/nSites),2),
      ")"
      ))

noquote(paste("Mean of B Scores: ",
              signif(mean(Bscores), 3),
              " (SE ",
              signif(sqrt(var(Bscores)/nSites),2),
              ")"
))

plot(as.numeric(substr(siteNames, 2, 6)), Ascores, type="l",
     xlab="Position in mitogenome",
     ylab="Score")
points(as.numeric(substr(siteNames, 2, 6)), Bscores, col="grey", type="l")
legend("top", lty=1, col=c("black","grey"), legend=c("Ascores","Bscores"))

##########################
# other ####
##########################
# mainDFfixed <- mainDF[mainDF$Position %in% fixedPos,]
# mainDFfixed$pop <- ifelse(mainDFfixed$Sample %in% pop0, "1", "2")
# head(mainDFfixed)

# 
# plotFixed <- function(x){
#   i <- fixedPos[x]
#   print(x)
#   datInt <- mainDFfixed[mainDFfixed$Position==i,]
#   datInt <- datInt[!is.na(datInt$ylog),]
#   datInt <- datInt[!is.na(datInt$xnqlogis),]
#        plot(ylog~xnqlogis,
#             xlim=c(0, 9.5),
#             ylim=c(-8,0),
#             main = i,
#             col = (datInt$pop=="1")+1,
#             pch = 19,
#             data=datInt)
# 
# 
# }
# plotFixed(2)
# plotFixed(80)
# 
# 
# par(mfrow=c(5,5))
# sapply(1:25, function(x) plotFixed(x))
# sapply(26:50, function(x) plotFixed(x))
# sapply(51:75, function(x) plotFixed(x))
# sapply(76:100, function(x) plotFixed(x))
# sapply(101:111, function(x) plotFixed(x))
# par(mfrow=c(1,1))
# 
# dim(mainDFfixed)
# head(mainDFfixed)
# # Pop 1 always has the lower allele freqs
# all(
#   sapply(unique(mainDFfixed$Position), function(x){
#   dff <- mainDFfixed[mainDFfixed$Position == x, ]
#   mean(dff[dff$pop=="1", 3], na.rm=T) < mean(dff[dff$pop=="2", 3], na.rm =T)
# })
# )
# 
# # flip allel frqs in pop 2
# mainDFfixed[mainDFfixed$pop == "2", 3] <- 1 - mainDFfixed[mainDFfixed$pop == "2", 3]
# # adjust qlogis vals
# mainDFfixed$ylog <- log(mainDFfixed$AltProp)
# mainDFfixed$ylog[is.infinite(mainDFfixed$ylog)] <- NA
# 
# head(mainDFfixed)
# models <- lapply(unique(mainDFfixed$Position), function(x){
#   dff <- mainDFfixed[mainDFfixed$Position == x, ]
#   plot(dff$ylog~dff$xnqlogis, main = x)
#   with(dff,{
#     lm(ylog~1+pop, offset=xnqlogis)
#   }
#   )
# })
# 
# summary(models[[3]])
# intercepts <- t(sapply(models, function(x) coef(summary(x))[,1]))
# interceptsAdj <- t(apply(intercepts, 1, function(x) c(x[1], sum(x))))
# head(intercepts)
# head(interceptsAdj)
# plot(intercepts)
# plot(interceptsAdj, asp=1, xlim=c(-12,-7), ylim=c(-12,-7))
# grid()
# abline(0,1)
# interceptsAdjHopper <- interceptsAdj # for use in parrot script
# 
# hist(apply(interceptsAdj, 1, max))
# meanEstFromFixed <- exp(apply(interceptsAdj, 1, max))
# shapiro.test(meanEstFromFixed) # not significantly different from normal
# # Mean estimate
# mean(meanEstFromFixed)*100
# # 95% confidence interval
# quantile(meanEstFromFixed, c(0.025, 0.975))*100
# 
# hist(meanEstFromFixed*100, main="Grasshopper (fixed SNPs)", xlab="Estimate of NUMT proportion (in percent)")
# abline(v=mean(meanEstFromFixed*100), lwd=2)
# abline(v=quantile(meanEstFromFixed*100, c(0.025, 0.975)), lty=2, lwd=2)
# legend("topleft", lty=c(1, 2), lwd=2, legend=c("Mean","95% CI"))
