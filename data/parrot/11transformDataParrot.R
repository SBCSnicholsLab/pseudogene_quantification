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
# 
# # population divergence ####
# 
# 
# dim(gtsHC)
# head(gtsHC)
# pc01 <- prcomp(t(gtsHC[,1:nSamp]))
# 
# # PC1 accounts for >89% of the variance, PC2 for <3%
# summary(pc01)
# str(pc01)
# pc01$x[,1:2]
# 
# # PCA separates two clusters of individuals
# plot(pc01$x[,1:2], asp = 3/90, main="PCA mitotypes")
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
# fixedIndAnyCov <- apply(gtsHC, 1, function(x) is.fixed.diff(x, pop0ind, pop1ind))
# sum(fixedIndAnyCov)
# 
# fixedPos <- paste0("S", sprintf("%05d", gtsHC$POS[fixedIndAnyCov]))
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
#   plot(ylog~xnqlogis,
#        xlim=c(0, 9.5),
#        ylim=c(-8,0),
#        main = i,
#        col = (datInt$pop=="1")+1,
#        pch = 19,
#        data=datInt)
#   
#   
# }
# plotFixed(2)
# plotFixed(80)
# 
# 
# 
# dim(mainDFfixed)
# head(mainDFfixed)
# # Pop 2 always has the lower allele freqs
# all(
#   sapply(unique(mainDFfixed$Position), function(x){
#     dff <- mainDFfixed[mainDFfixed$Position == x, ]
#     mean(dff[dff$pop=="1", 3], na.rm=T) > mean(dff[dff$pop=="2", 3], na.rm =T)
#   })
# )
# 
# # flip allel frqs in pop 1
# mainDFfixed[mainDFfixed$pop == "1", 3] <- 1 - mainDFfixed[mainDFfixed$pop == "1", 3]
# # adjust qlogis vals
# mainDFfixed$ylog <- log(mainDFfixed$AltProp)
# 
# mainDFfixed <- mainDFfixed[!is.infinite(mainDFfixed$ylog),]
# 
# head(mainDFfixed)
# models <- lapply(unique(mainDFfixed$Position), function(x){
#   dff <- mainDFfixed[mainDFfixed$Position == x, ]
#   plot(dff$ylog~dff$xnqlogis, main = x)
#   if(length(table(dff$pop))!=2) return(NULL)
# 
#   with(dff,{
#     lm(ylog~1+pop, offset=xnqlogis)
#   }
#   )
# })
# which(sapply(models, is.null))
# models <- models[-which(sapply(models, is.null))]
# summary(models[[3]])
# intercepts <- t(sapply(models, function(x) coef(summary(x))[,1]))
# interceptsAdj <- t(apply(intercepts, 1, function(x) c(x[1], sum(x))))
# head(intercepts)
# head(interceptsAdj)
# plot(intercepts)
# plot(interceptsAdj, asp=1)
# grid()
# abline(0,1)
# 
# # whichever intwrcept is the higher one, is the estimate!
# estDiverg <- apply(interceptsAdj, 1, max)
# hist(exp(estDiverg)*100)
# 
# 
# 
# 
# hist(diffs)
# 
# coef(summary(models[[1]]))[3,1]
