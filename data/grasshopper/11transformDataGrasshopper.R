######################################
# Prepare data for general method ####
######################################

################
# mapping depths

setwd("~/git_repos/pseudogene_quantification/data/grasshopper/")
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

# 
# # population divergence ####
# 
# 
# dim(gtsHC)
# head(gtsHC)
# pc01 <- prcomp(t(gtsHC[,1:nSamp]))
# 
# # PC1 accounts for >60% of the variance, PC2 for <3%
# summary(pc01)
# str(pc01)
# pc01$x[,1:2]
# 
# # PCA separates two clusters of individuals
# plot(pc01$x[,1:2], asp = 6/64, main="PCA mitotypes")
# 
# # Which sites have fixed differences between the groups? ####
# 
# # Get sample names of PCA rgoups
# pop0 <- rownames(pc01$x)[pc01$x[,1]>0]
# pop1 <- rownames(pc01$x)[pc01$x[,1]<0]
# 
# # Get indices
# pop0ind <- which(names(gtsHC) %in% pop0)
# pop1ind <- which(names(gtsHC) %in% pop1)
# 
# # A function to check whether there is sharing of alleles
# #  (for an individual locus)
# is.fixed.diff <- function(x, ind0, ind1){
#   !any(x[ind0] %in% x[ind1])
# }
# 
# # The first locus in gtsHC is not fixed:
# is.fixed.diff(gtsHC[1,], pop1ind, pop0ind)
# 
# cov15 <- c(names(which(mappingProp > 0.002)), "POS")
# 
# 
# names(gtsHC) %in% cov15
# gtsCov15 <- gtsHC[, names(gtsHC) %in% cov15]
# 
# pop0.15 <- pop0[pop0 %in% names(gtsCov15)]
# pop1.15 <- pop1[pop1 %in% names(gtsCov15)]
# pop0ind.15 <- which(names(gtsCov15) %in% pop0.15)
# pop1ind.15 <- which(names(gtsCov15) %in% pop1.15)
# 
# 
# head(gtsCov15)
# fixedInd <- apply(gtsCov15, 1, function(x) is.fixed.diff(x, pop0ind.15, pop1ind.15))
# sum(fixedInd)
# 
# fixedPos <- paste0("S", sprintf("%05d", gtsHC$POS[fixedInd]))
# 
# 
# mainDFfixed <- mainDF[mainDF$Position %in% fixedPos,]
# mainDFfixed$pop <- ifelse(mainDFfixed$Sample %in% pop0, "1", "2")
# head(mainDFfixed)
# 
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
