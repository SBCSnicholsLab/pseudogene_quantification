######################################
# Prepare data for general method ####
######################################

################
# mapping depths

setwd("~/git_repos/pseudogene_quantification/data/parrot/")
dir()
nSamp = 22
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
hist(mappingProp, breaks=100)


#################
# Gentotype calls

gts <- read.table("genotypes.csv")
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
# 7784*4 = 31136

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
# long(er) DF with GT calls
gtDF <- data.frame(sample = sampleIndex, pos = rep(gtsHC$POS, nSamp),
                   gt = unlist(gtsHC[,1:nSamp]) + 3)
head(gtDF)
str(gtDF)

# Depth of the alt allels (summed)
siteDeps <- rowSums(allDF[,3:4])
altRatio <- allDF$alt / siteDeps



# DF that wih all releavant data for the anaysis:
mainDF <- data.frame(Sample = as.character(gtDF$sample),
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

# with(mainDF,{
#   plot(ylog~xnqlogis, pch = '.',
#        main="Podi data",
#        xlab="logit of mapping proportion",
#        ylab="log rel minor allele freq")
# })

# Write out data
#write.table(mainDF, "transformedData.csv")



#################################################
# Diverged populations and fixed differences ####
#################################################
# 

# population divergence ####

# samples and SRA ids from NCBI
idTab <- as.data.frame(rbind(
c("SRR6214420", "B41207"),
c("SRR6214421", "B4471"),
c("SRR6214422", "B2117"),
c("SRR6214423", "B22946"),
c("SRR6214424", "B18697"),
c("SRR6214425", "B20203"),
c("SRR6214426", "B2966"),
c("SRR6214427", "B32871"),
c("SRR6214428", "B22948"),
c("SRR6214429", "B28284"),
c("SRR6214430", "R7511"),
c("SRR6214431", "HLW99"),
c("SRR6214432", "HLW77"),
c("SRR6214433", "HLW90"),
c("SRR6214434", "B49688"),
c("SRR6214435", "B49747"),
c("SRR6214436", "B46433"),
c("SRR6214437", "B47715"),
c("SRR6214438", "B9297"),
c("SRR6214439", "HLW7303"),
c("SRR6214440", "B51936"),
c("SRR6214441", "B54105")
))
names(idTab) <- c("SRA", "ID")

dim(gtsHC)
head(gtsHC)
pc01 <- prcomp(t(gtsHC[,1:nSamp]))

# PC1 accounts for >89% of the variance, PC2 for <3%
summary(pc01)
str(pc01)
pc01$x[,1:2]
pcCoords <- data.frame(SRA=rownames(pc01$x[,1:2]), pc01$x[,1:2])
pcCoords <- merge(idTab, pcCoords, by="SRA")
# PCA separates two clusters of individuals
plot(pcCoords[,3:4], main="PCA mitotypes")
plot(pcCoords[,3:4], main="PCA mitotypes", asp=3/90)

# Same groups as in McElroy 2018:
as.character(pcCoords[pcCoords[,3]<0,]$ID)
as.character(pcCoords[pcCoords[,3]>0,]$ID)


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

fixedIndAnyCov <- apply(gtsHC, 1, function(x) is.fixed.diff(x, pop0ind, pop1ind))
sum(fixedIndAnyCov)

fixedPos <- paste0("S", sprintf("%05d", gtsHC$POS[fixedIndAnyCov]))


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
head(allCountsFixedN)
allCountsFixedN <- allCountsFixedN[,-c(7:10)]
# write.table(allCountsFixedN, "parrotFixed.csv")
# # update the major allele matrix
# oneTrue <- function(x){
#   a = c(F, F, F, F)
#   a[x] <- T
#   a
# }
# # oneTrue(2)
# # This is only 6k rows, so we can use sapply.
# majMat <- as.matrix(t(sapply(1:nrow(allCountsFixedN), function(x) {
#   c(oneTrue(allCountsFixedN[x, ]$A),
#     oneTrue(allCountsFixedN[x, ]$B)
#   )})))
# # only one major allele now?
# all(rowSums(majMat) == 2)
# # yes, there are two because there's one major allele in A and one in B
# colnames(majMat) <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4")
# dim(majMat)
# dim(allCountsFixedN)
# head(allCountsFixedN)
# allCountsFixedN <- data.frame(allCountsFixedN[, -c(7:10)], majMat)
# head(allCountsFixedN)
# allCountsFixedN$allDep <- rowSums(allCountsFixedN[,3:6])
# head(allCountsFixedN)
# # get numbers for each locus and allele aMaj and bAlt are sometimes indentical,
# # but only if there were only two alleles.
# aMaj <- rowSums(allCountsFixedN[,3:6] * allCountsFixedN[,12:15])
# bMaj <- rowSums(allCountsFixedN[,3:6] * allCountsFixedN[,16:19])
# aAlt <- rowSums(allCountsFixedN[,3:6] * !allCountsFixedN[,12:15])
# bAlt <- rowSums(allCountsFixedN[,3:6] * !allCountsFixedN[,16:19])
# allCountsFixedN <- data.frame(allCountsFixedN, aMaj, aAlt, bMaj, bAlt)
# head(allCountsFixedN)


# Richard's scores ###################################################
######################################################################

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

ABdf <- data.frame(sA = Ascores, sB=Bscores)
plot(sB~sA, data=ABdf)
plot(sB~sA, data=ABdf[-161,]) # outlier at fixed site 161 removed

lm0 <- lm(sB~sA, data=ABdf)
lm1 <- lm(sB~sA, data=ABdf, offset = sA)
lm2 <- lm(sB~sA, data=ABdf[-161,])
lm3 <- lm(sB~sA, data=ABdf[-161,], offset = sA)

summary(lm0) # all fixed sites
summary(lm1) # all fixed sites with offset
summary(lm2) # -161
summary(lm3) # -161 with offset

abline(lm0)
abline(lm2) # more like what we expect

abline(0, 1, lty = 2)
abline(2.5e-6, 0.9, lty = 2)
par(mfrow=c(2,2))
plot(lm0)
plot(lm2) # 161 removed, still not great

par(mfrow=c(1,1))


mappingDepths <- data.frame(nMapped, pop=rep("A", length(nMapped)), stringsAsFactors = F)
mappingDepths$pop[pop0ind]
mappingDepths$pop[pop1ind] <- "B"
summary(lm(nMapped ~ pop, data=mappingDepths))
#plot(lm(nMapped ~ pop, data=mappingDepths))



noquote(paste("Mean of A Scores: ",
              signif(mean(Ascores[-161]), 3),
              " (SE ",
              signif(sqrt(var(Ascores[-161])/(nSites-1)),2),
              ")"
))

noquote(paste("Mean of B Scores: ",
              signif(mean(Bscores[-161]), 3),
              " (SE ",
              signif(sqrt(var(Bscores)/(nSites-1)),2),
              ")"
))
hist(Ascores[-161], col = "#0000FF40", main = "Parrot", xlab="Scores", ylim=c(0, 90), xlim=c(0, 1/1000))
abline(v=mean(Ascores[-161]), lwd=2)
abline(v=mean(Ascores[-161]) + c(-1, 1) * sqrt(var(Ascores[-161])/(nSites-1)), lty=2)
# intervals overlap
# hist(Bscores[-161], col = "#FF000040", add=T)
# abline(v=mean(Bscores[-161]), lwd=2)
# abline(v=mean(Bscores[-161]) + c(-1, 1) * sqrt(var(Bscores[-161])/(nSites-1)), lty=2)
# legend("topright", fill=c("#0000FF40","#FF000040"), legend=c("A","B"))

