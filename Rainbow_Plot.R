################################################################################
# A method to estimate the nuclear frequency of exogenous pseudogenes from     #
#   population-level low-pass sequencing data                                  #
################################################################################

library(lme4) 




#' Rainbow Plot
#' Generate a plot to estimate the proportion of the nuclear genome
#' made up of a particular vagrant sequence. The function produces the plots
#' intercept estimate, and mapping depth estimate as described by
#' Becher & Nichols (2021) citation.
#' 
#' @param data 
#' @param AltProp 
#' @param Position 
#' @param ylog 
#' @param xnqlogis 
#' @param Sample 
#' @param nloci 
#' @param minWt 
#' @param maxFreq 
#' @param minSamples 
#' @param filterHard 
#' @param seed 
#' @param main 
#' @param printout 
#'
#' @return
#' @export
#'
#' @examples
rainbowPlot <- function(data,
                        AltProp = data$AltProp,
                        Position = factor(data$Position),
                        ylog = data$ylog,
                        xnqlogis = data$xnqlogis,
                        Sample = factor(data$Sample),
                        nloci = 200,
                        minWt = 0.01,
                        maxFreq = 0.8,
                        minSamples = 10,
                        filterHard = TRUE,
                        seed,
                        main="",
                        printout = TRUE
){
  
  if( any( c( !is.vector(AltProp), 
              !is.factor(Position),
              !is.vector(xnqlogis),
              !is.vector(ylog),
              !is.factor(Sample)
  )
  )
  ) stop("You must either supply a dataframe containing a factors Position & Sample (or vectors which can be coerced to a factor), \n 
              \t plus vectors AltProp, ylog & xnqlogis \n 
              \t or specify them as separate Arguments")
  
  if( !require(lme4)
  ) stop("Install package lme4 before running this function")
  
  # Save the function call
  funcCall <- sys.call()
  
  # reconstruct data.frame having coerced Sample and Position to factors
  data <- data.frame(Sample, Position, AltProp, ylog, xnqlogis)
  
  # Remove the rows of data.frame with NAs 
  goodDat <- subset(data,complete.cases(data))
  
  # Remove abnormally high values
  goodDat <- subset(goodDat, AltProp < maxFreq)
  
  
  # Find how many individuals each SNP has been scored in
  tt <- table(goodDat$Position)
  # get the names of SNPs occurring in more than minSamples
  goodNames <- names(tt)[tt>minSamples]
  
  # Exclude loci found in too few samples
  goodDat <- subset(goodDat, Position %in% goodNames)
  
  
  # set random number seed if requested (ie. if seed is specified in the function call)
  if(hasArg(seed)) set.seed(seed)
  
  # Find the average allele frequency at each site
  wweights <- tapply(goodDat$AltProp, factor(goodDat$Position), mean)
  # exclude loci with average frequency less than minWt
  wweights[wweights < minWt] <- 0
  
  # sample a set of high frequency loci (sampling proportional to the freq)
  loci <- sample(names(wweights), nloci, prob = wweights)
  
  
  # Calculate the raw slopes for each locus (to identify any outlying values)
  slopes <- rep(0, nloci)
  for (i in 1:nloci) slopes[i] <- coef(lm(ylog ~ xnqlogis, 
                                          data = subset(goodDat, Position == loci[i])
  )
  )[2]
  
  # Use the boxplot.stats function to identifying SNPs with outlying slopes
  # Hard filtering excludes upper and lower quartiles
  # otherwise the boxplot default criterion is used (familiar as the end of whiskers)
  bpsSlopes <- boxplot.stats(slopes)$stats
  if (filterHard) {
    rogueSNPs <- which(slopes > bpsSlopes[4] | slopes < bpsSlopes[2])} else {
      rogueSNPs <- which(slopes > bpsSlopes[5] | slopes < bpsSlopes[1])}
  goodLoci <- loci[-rogueSNPs]
  
  # remove the rogue loci from the data.frame
  goodDat <- subset(goodDat,Position %in% goodLoci)
  
  
  
  # fit a 1:1 line plus intercept plus sample random effect  
  lmod5 <- lmer(ylog ~ 0 + Position + (1 | Sample),
                offset = xnqlogis,
                data = goodDat)
  
  
  intercepts5 <- summary(lmod5)$coefficients[,1]
  
  SEs5 <- summary(lmod5)$coefficients[,2]
  
  
  
  # rank intercepts (ranking is used to select SNPs colour on the plot & associated lines)
  iranks <- rank(intercepts5, ties.method = "random")
  # create a data.frame associating the name of the SNP with its rank
  rankDF <- data.frame(Position=substr(names(iranks), 9, 20), rank=iranks)
  # use the merge function to add these values the appropriate row in goodDat
  goodDat <- merge(goodDat, rankDF, by="Position")
  
  # Put the raw data points for selected SNPs onto the rainbow plot  
  maxx <- max(c(10, goodDat$xnqlogis))
  miny <- min(c(-10,max(intercepts5)))
  plot(ylog~xnqlogis, col=rainbow(max(rank)*1.4)[rank], 
       pch=1, xlim=c(0,maxx), ylim=c(miny,0),
       xlab="(Un-mapped data / mapped)",
       ylab="Allele frequency",
       main=main,
       xaxt = 'n',
       yaxt = 'n',
       data = goodDat
  )
  
  # add the axes 
  axis(1, at = log(2^(0:7*2)), labels = 2^(0:7*2))
  abline(v = log(2^(0:7*2)), col = "grey", lty=3)
  ll <- expression("1", "10"^-1, "10"^-2, "10"^-3, "10"^-4, "10"^-5, "10"^-6)
  axis(2, at = log(10^(0:-6)), labels = ll)
  abline(h = log(10^(0:-6)), col = "grey", lty=3)
  
  # Add key lines
  abline(h=max(intercepts5))
  abline(v=max(goodDat$xnqlogis))
  abline(max(intercepts5), 1,  lty=2)
  
  # Get intercepts and standard errors
  estIntlog <- max(intercepts5)
  estIntSE <- SEs5[which(intercepts5 == max(intercepts5))][1] # S.E. in log space
  
  # Convert CI to real values
  estInt <- plogis(unname(c(intEst=estIntlog,
                            intLo=estIntlog - (1.96 * estIntSE),
                            intUp=estIntlog + (1.96 * estIntSE)
  ))
  )
  
  # Add annotation to the plot
  estDep <- plogis(-max(goodDat$xnqlogis))
  text(c(2, 7), c(0,0), c(
    paste0("intercept est: ", round(estInt[1]*100, digits = 3), "%\n",
           "(", round(estInt[2]*100, digits = 3), "%-",
           round(estInt[3]*100, digits = 3), "%)"),
    paste0("mapping depth est: ", round(estDep*100, digits = 3), "%")
  )
  
  )
  
  # Add the trend lines to the rainbow plot
  sapply(1:length(intercepts5), function(x){
    abline(intercepts5[x], 1, col = rainbow(max(iranks)*1.4, alpha = 0.2)[iranks][x])
  })
  
  
  numgoodloci <- length(goodLoci)
  
  if (printout) {
    cat('Intercept based on ', numgoodloci, 'SNP loci \n')
    cat('Estimate: ', signif(estInt[1],3), '\n')
    cat('Confindence Interval: ',
        signif(estInt[2],3),
        '-',
        signif(estInt[3],3),
        '\n'
    )
    cat('Mapping depth estimate: ',signif(estDep,3), '\n')
    cat('Function call \n', deparse(funcCall))
  }
  
  # Return the estimated values as an (initially) invisible object for further use
  invisible(list(intercepts = c(intercept.est=estInt[1],
                                intercept.est.lo=estInt[2],
                                intercept.est.up=estInt[3]),
                 depth.est = estDep,
                 num.loci = numgoodloci,
                 func.params = funcCall)
  )
  
}


# Get data online ####
# human
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/human/transformedData.csv",
              destfile = "human.csv")
humanDF <- read.table("human.csv", stringsAsFactors = F)

#hopper
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/grasshopper/transformedData.csv",
              destfile = "hopper.csv")
hopperDF <- read.table("hopper.csv", stringsAsFactors = F)

#parrot
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/parrot/transformedData.csv",
              destfile = "parrot.csv")
parrotDF <- read.table("parrot.csv", stringsAsFactors = F)

# Alternatively, read in data from repo ####
# setwd("~/git_repos/pseudogene_quantification/")
# # human data set
# humanDF <- read.table("data/human/transformedData.csv", stringsAsFactors = F)
# # parrot data set
# parrotDF <- read.table("data/parrot/transformedData.csv", stringsAsFactors = F)
# # grasshopper data set
# hopperDF <- read.table("data/grasshopper/transformedData.csv", stringsAsFactors = F)


rainbowPlot(humanDF, seed = 12345, main="Human")

rainbowPlot(hopperDF, seed = 12345, main="Grasshopper")

rainbowPlot(parrotDF, seed = 12345, main="Parrot")


