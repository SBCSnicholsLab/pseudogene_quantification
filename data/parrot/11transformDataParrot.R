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
# 
# # population divergence ####
# 
# # samples and SRA ids from NCBI
# idTab <- as.data.frame(rbind(
# c("SRR6214420", "B41207"),
# c("SRR6214421", "B4471"),
# c("SRR6214422", "B2117"),
# c("SRR6214423", "B22946"),
# c("SRR6214424", "B18697"),
# c("SRR6214425", "B20203"),
# c("SRR6214426", "B2966"),
# c("SRR6214427", "B32871"),
# c("SRR6214428", "B22948"),
# c("SRR6214429", "B28284"),
# c("SRR6214430", "R7511"),
# c("SRR6214431", "HLW99"),
# c("SRR6214432", "HLW77"),
# c("SRR6214433", "HLW90"),
# c("SRR6214434", "B49688"),
# c("SRR6214435", "B49747"),
# c("SRR6214436", "B46433"),
# c("SRR6214437", "B47715"),
# c("SRR6214438", "B9297"),
# c("SRR6214439", "HLW7303"),
# c("SRR6214440", "B51936"),
# c("SRR6214441", "B54105")
# ))
# names(idTab) <- c("SRA", "ID")
# 
# dim(gtsHC)
# head(gtsHC)
# pc01 <- prcomp(t(gtsHC[,1:nSamp]))
# 
# # PC1 accounts for >89% of the variance, PC2 for <3%
# summary(pc01)
# str(pc01)
# pc01$x[,1:2]
# pcCoords <- data.frame(SRA=rownames(pc01$x[,1:2]), pc01$x[,1:2])
# pcCoords <- merge(idTab, pcCoords, by="SRA")
# # PCA separates two clusters of individuals
# plot(pcCoords[,3:4], main="PCA mitotypes")
# plot(pcCoords[,3:4], main="PCA mitotypes", asp=3/90)
# 
# # Same groups as in McElroy 2018:
# as.character(pcCoords[pcCoords[,3]<0,]$ID)
# as.character(pcCoords[pcCoords[,3]>0,]$ID)
# 
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
# plotFixed(3)
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
# # whichever intercept is the higher one, is the estimate!
# estDiverg <- apply(interceptsAdj, 1, max)
# hist(exp(estDiverg)*100)
# 
# 
# 
# 
# # plot comparing to grasshopper values
# plot(interceptsAdjHopper, asp=1, xlim=c(-10,-7), ylim=c(-13,-7),
#      xlab= "Relative allele frequency in population A",
#      ylab= "Relative allele frequency in population B",
#      xaxt = 'n',
#      yaxt = 'n',
#      type='n'
# )
# 
# axis(1, at = log(c(10^-5, 3*10^-5, 10^-4, 3*10^-4, 10^-3)),
#      labels = expression("10"^-5, "3×10"^-5, "10"^-4, "3×10"^-4, "10"^-3))
# abline(v = log(c(10^-5, 3*10^-5, 10^-4, 3*10^-4, 10^-3)), col = "grey", lty=3)
# 
# axis(2, at = log(c(3*10^-6, 10^-5, 3*10^-5, 10^-4, 3*10^-4, 10^-3)),
#      labels = expression("3×10"^-6, "10"^-5, "3×10"^-5, "10"^-4, "3×10"^-4, "10"^-3))
# abline(h = log(c(3*10^-6, 10^-5, 3*10^-5, 10^-4, 3*10^-4, 10^-3)), col = "grey", lty=3)
# 
# 
# 
# points(interceptsAdjHopper, col = 1)
# points(interceptsAdj, col = 2)
# 
# abline(0,1, lty=2)
# legend("bottomleft", pch=1, col = c(1,2), legend = c("Grasshopper", "Parrot"))
