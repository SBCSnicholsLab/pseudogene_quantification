######################################
# Prepare data for general method ####
######################################

################
# mapping depths

# adjust directory
setwd("~/git_repos/pseudogene_quantification/data/parrot/")
nSamp <- 22

# a utility function to rename the mapping depth vals
getNums <- function(x){
  a <- x[,3]
  names(a) <- sapply(x[,1], function(y){
    strsplit(y, "j")[[1]][1]
  })
  return(a)
}
nMapped <- read.table("bpMapped", sep = "\t", stringsAsFactors = F)
#head(nMapped)

nTot <- read.table("bpTotal", sep = "\t", stringsAsFactors = F)
#head(nTot)

nMapped <- getNums(nMapped)
nTot <- getNums(nTot)

#################
# Gentotype calls

gts <- read.table("genotypes.csv",header = T, check.names = F)
#head(gts)

# remove lowest coverage samples
gtsHC <- gts # none removed

# visualise missingness
#image(as.matrix(gtsHC[,1:nSamp]) == -1)

# remove individuals with missing data
# indMissing <- apply(gtsHC, 2, function(x) -1 %in% x)
# 
# 
# gtsHC <- gtsHC[,!indMissing]
#head(gtsHC)

# # remove any site with missing data
posMissing <- gtsHC$POS[apply(gtsHC, 1, function(x) -1 %in% x)]
length(posMissing) # only 67 sites. Remove all
gtsHC <- gtsHC[!(gtsHC$POS %in% posMissing),]

nindFilt <- ncol(gtsHC)-1
###############
# Allele counts

allCounts <- read.table("alleleCounts.csv")
#head(allCounts)
allCountsHC <- allCounts # opportunity to remove individuals if desired

# remove sites with missing alleles identified above
allCountsHC <- allCountsHC[!(allCountsHC$position %in% posMissing),]
#head(allCountsHC)

# split allele counts into separate DFs, one per allele
gtList <- lapply(1:4, function(x) {
  allCountsHC[allCountsHC$allele==x,]
})


# flatten into four vectors
#rm(c) # In case somebody named a variable c, which would shadow the c function.
gtVectors <- lapply(gtList, function(x){
  do.call(c,x[,1:nindFilt])
})
#index of samples
sampleIndex <- rep(colnames(gtsHC)[1:nindFilt], each=nrow(gtList[[1]]))


# long(er) data frame for allele counts with sample index
allDF <- data.frame(sample=sampleIndex,
                    pos = rep(gtsHC$POS, nindFilt),
                    ref=gtVectors[[1]],
                    alt=gtVectors[[2]] + gtVectors[[3]] + gtVectors[[4]])
#head(allDF)
#str(allDF)

# Depth of the alt allels (summed)
siteDeps <- rowSums(allDF[,3:4])
altRatio <- allDF$alt / siteDeps



# DF that wih all releavant data for the anaysis:
mainDF <- data.frame(Sample = as.character(allDF$sample),
                     Position=allDF$pos,
                     AltProp=altRatio,
                     DP=allDF$ref+allDF$alt,
                     stringsAsFactors = F, row.names = NULL)
#head(mainDF)
# Add mapping rate, make a DF so we can use "merge"
#mappingProp
propDF <- data.frame(Sample = names(nMapped),
                     nMapped,
                     nTot,
                     stringsAsFactors = F)
#head(propDF)
#head(mainDF)
mainDF <- merge(mainDF, propDF, by="Sample")

#head(mainDF)
mainDF$ylog <- log(mainDF$AltProp)
mainDF$ylog[is.infinite(mainDF$ylog)] <- NA # replace -Inf by NA, makes plotting easier




#head(mainDF)

# remove NA lines
mainDF <- mainDF[!is.na(mainDF$AltProp),]
mainDF <- mainDF[!is.na(mainDF$ylog),]


# Write out data
write.table(mainDF, "transformedData.csv")




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

plot(pcCoords[,3:4], main="PCA mitotypes", asp=3/90)

# Same groups as in McElroy et al. 2018:
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

#library(vagrantDNA)
#divEst(allCountsFixedN)
