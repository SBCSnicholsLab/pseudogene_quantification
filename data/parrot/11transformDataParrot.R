


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

# plit allele counts into separate DFs, one per allele
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
                    A1=gtVectors[[1]], A2=gtVectors[[2]], A3=gtVectors[[3]],
                    A4=gtVectors[[4]])
head(allDF)
str(allDF)
# long(er) DF with GT calls
gtDF <- data.frame(sample = sampleIndex, pos = rep(gtsHC$POS, nSamp),
                   gt = unlist(gtsHC[,1:nSamp]) + 3)
head(gtDF)
str(gtDF)


# depth of the called allele. This does a lot of indexing,
#  picking the counts of the allele that was called.
callDeps <- sapply(1:nrow(gtDF), function(x) {
  allDF[x, gtDF[x, 3]]
})


# Depth of the alt allels (summed)
siteDeps <- rowSums(allDF[,3:6])
altDeps <- siteDeps - callDeps
altRatio <- altDeps / siteDeps

siteSums <- tapply(siteDeps, allDF$pos, sum)
altSums <- tapply(altDeps, allDF$pos, sum)

# Allele frequency spectrum
sumRatios <- as.vector(altSums/siteSums)
hist(sumRatios, breaks=200)
abline(v=0.01)

rareVariants <- names(sumRatios[sumRatios < 0.01])
nrow(gtsHC)


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











# selection of loci ####

# overall mapping depths ####
# numt alleles are low, what are the overall mapping depths?
head(allCountsHC)
siteSUMS <- tapply(1:nrow(allCountsHC), allCountsHC$position, function(x){
  colSums(allCountsHC[x, 1:nSamp])
})
siteSUMS <- as.matrix(do.call(rbind, siteSUMS))
dim(siteSUMS)
head(siteSUMS)
barplot(colMeans(siteSUMS), las=3, log="y")
abline(h=c(10, 50, 100, 200, 500))

# population divergence ####


dim(gtsHC)
head(gtsHC)
pc01 <- prcomp(t(gtsHC[,1:nSamp]))
summary(pc01)
str(pc01)
pc01$x[,1:2]
plot(pc01$x[,1:2], asp = 2/90)
plot(pc01$x[,1:2])
abline(h=c(0.5, -1))
which(pc01$x[,2] < -1)
which(pc01$x[,2] > 0.5)
which(pc01$x[,1] > 0)

plot(pc01$rotation[,1])
which(pc01$rotation[,1] > 0.15)
hist(mappingProp, breaks=20)
pop0 <- rownames(pc01$x)[pc01$x[,1]>0]
pop1 <- rownames(pc01$x)[pc01$x[,1]<0]

pop0ind <- which(names(gtsHC) %in% pop0)
pop1ind <- which(names(gtsHC) %in% pop1)

is.fixed.diff <- function(x, ind0, ind1){
  !any(x[ind0] %in% x[ind1])
}
is.fixed.diff(gtsHC[1,], pop1ind, pop0ind)

fixedInd <- apply(gtsHC, 1, function(x) is.fixed.diff(x, pop0ind, pop1ind))
sum(fixedInd)

fixedPos <- paste0("S", sprintf("%05d", gtsHC$POS[fixedInd]))


mainDFfixed <- mainDF[mainDF$Position %in% fixedPos,]
mainDFfixed$parrotpop <- mainDFfixed$Sample %in% pop0

a <- numeric(0)
for( j in 1:length(fixedPos)) {
  i <- fixedPos[j]
  if(i %in% mainDFfixed$Position){
    with(mainDFfixed[mainDFfixed$Position==i,],
         plot(ylog~xnqlogis,
              xlim=c(0, 7.5),
              ylim=c(-8,0),
              main = i,
              col = (parrotpop==T)+1,
              pch = 19))
    mod00 <- lm(ylog~0+parrotpop, data = mainDFfixed[mainDFfixed$Position==i,], offset=xnqlogis)
    #print(summary(mod00))
    abline(coef(mod00)[1], 1)
    abline(coef(mod00)[2], 1)
    phat <- plogis(diff(coef(mod00)))
    logc <- mean(coefficients(mod00)-c(log(1-phat),log(phat)))
    a[j] <- (exp(logc))
  }
  else
    print(i, " is missing")
}


plotParrotFixed <- function(x){
  i <- fixedPos[x]
  print(x)
  datInt <- mainDFfixed[mainDFfixed$Position==i,]
  datInt <- datInt[!is.na(datInt$ylog),]
  datInt <- datInt[!is.na(datInt$xnqlogis),]
       plot(ylog~xnqlogis,
            xlim=c(0, 9.5),
            ylim=c(-8,0),
            main = i,
            col = (datInt$parrotpop==T)+1,
            pch = 19,
            data=datInt)
  #print(datInt)
  if( nrow(datInt) > 0){
    mod00 <- lm(ylog~0+parrotpop, data = datInt, offset=xnqlogis)
    #print(summary(mod00))
    if(all(is.finite(coef(mod00)))){
      abline(coef(mod00)[1], 1)
      abline(coef(mod00)[2], 1)
      phat <- plogis(diff(coef(mod00)))
      logc <- mean(coefficients(mod00)-c(log(1-phat),log(phat)))
      return(exp(logc)*100)
    } else {
      return(NA)
    }
  } else return(NA)
  
}
plotParrotFixed(1)
plotParrotFixed(2)
a <- sapply(1:179, plotParrotFixed)
sum(is.na(a))
hist(a, breaks = 50)
