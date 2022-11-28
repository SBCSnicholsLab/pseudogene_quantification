######################################
# Prepare data for general method ####
######################################

################
# mapping depths

# adjust directory
setwd("~/git_repos/pseudogene_quantification/data/human02/")


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
nSamp <- ncol(gts)-1
# remove lowest coverage samples
gtsHC <- gts # none removed

# visualise missingness
#image(as.matrix(gtsHC[,1:nSamp]) == -1)
# not so much missing

# remove individuals with missing data
#indMissing <- apply(gtsHC, 2, function(x) -1 %in% x)
#gtsHC <- gtsHC[,!indMissing]
#head(gtsHC)

# # remove any site with missing data
# posMissing <- gtsHC$POS[apply(gtsHC, 1, function(x) -1 %in% x)]
# gtsHC <- gtsHC[,!gtsHC$POS %in% posMissing]

nindFilt <- ncol(gtsHC)-1
###############
# Allele counts

allCounts <- read.table("alleleCounts.csv")
#head(allCounts)
allCountsHC <- allCounts # opportunity to remove individuals if desired

# remove sites with missing alleles identified above
#allCountsHC <- allCountsHC[,!c(indMissing, F)]
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
#str(gtVectors)

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
head(mainDF)

