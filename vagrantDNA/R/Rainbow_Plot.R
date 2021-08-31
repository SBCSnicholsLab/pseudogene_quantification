################################################################################
# A method to estimate the nuclear frequency of exogenous pseudogenes from     #
#   population-level low-pass sequencing data                                  #
################################################################################

library(lme4)




#' Rainbow Plot.
#'
#' Generate a plot to estimate the proportion of the nuclear genome
#' made up of a particular vagrant sequence. The function produces the plots
#' intercept estimate, and mapping depth estimate described by
#' Becher & Nichols (2021) citation.
#'
#' @param data A data.frame with at least the following columns.
#' \describe{
#'   \item{AltProp}{A numberic vector giving the proportion of reads mapping to the exogenous genome}
#'   \item{Position}{A factor (or structure that can be coerced to a factor),
#'   giving a unique name for each SNP location.}
#'   \item{Sample}{A factor (or structure that can be coerced to a factor),
#'   giving a unique name for each sample.}
#'   \item{ylog}{A numeric vector giving the log(p), where p is the relative frequency of the SNP allele(s) not found in
#'   the exogenous genome.}
#'   \item{xnqlogis}{A numeric vector giving log(m/N);
#'   where m is the number of reads mapping to the exogenous genome and
#'   N is the remaining reads.}
#'   }
#' @param nloci The number of loci to be selected for the analysis. Default, 400.
#' @param minWt The minimum average allele frequency of SNPs to be included in the analysis.
#' Default, 0.01.
#' @param maxFreq The maximum allele frequency of an individual observation to be included
#' in the analysis. Default, 0.7.
#' @param minSamples The minimum number of samples in which a SNP should be called in order to
#' be included in the analysis.
#' @param filterHard If filterHard is TRUE, the SNP loci with slopes in the outer quartiles are
#' discarded. Otherwise the outliers identified by the default method of boxplot.stats function
#' are discarded. Default, TRUE.
#' @param seed Random number seed.
#' @param title User-supplied title for the rainbow plot.
#' @param printout If printout is TRUE, the function prints the estimates. Default, TRUE.
#' @return An invisible list with the following elements.
#' \describe{ 
#'  \item{$intercepts}{A vector giving the intercept estimate, 
#'  the lower and upper 95\% confidence interval bounds.}
#'  \item{$depth.est}{A crude upper limit; the minimum over samples 
#'  of the proportion of reads mapping to the exogenous genome.}
#'  \item{$num.loci}{The number of loci remaining after filtering,
#'  used to obtain the intercept estimates.}
#'  \item{$func.params}{A record of the function call}
#' }
#' @export
#'
#' @examples
#'
#'
rainbowPlot <- function(data,
                      nloci = 400,
                      minWt = 0.01,
                      maxFreq = 0.7,
                      minSamples = 10,
                      filterHard = TRUE,
                      seed,
                      title = "",
                      printout = TRUE
                      ){
  # force Position and Sample to become factors
  data$Position <- factor(data$Position)
  data$Sample <- factor(data$Sample)
  
  # check the dataframe has appropriate columns
  if( any( c( !is.vector(data$AltProp),
              !is.factor(data$Position),
              !is.vector(data$xnqlogis),
              !is.vector(data$ylog),
              !is.factor(data$Sample)
              )
           )
      ) stop("Supply a dataframe containing a factors Position & Sample (or vectors which can be coerced to a factor), \n
              \t plus vectors AltProp, ylog & xnqlogis")

  # check lme4 is installed
  if( !require(lme4)
      ) stop("Install package lme4 before running this function")

  # Save the function call
  funcCall <- sys.call()

  # Remove the rows of data with NAs
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

  # Find the average allele frequency at each Position
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
  goodDat <- subset(goodDat, Position %in% goodLoci)



  # fit a 1:1 line plus intercept plus sample random effect
  lmod5 <- lmer(ylog ~ 0 + Position + (1 | Sample),
                offset = xnqlogis,
                data = goodDat)


  intercepts5 <- summary(lmod5)$coefficients[,1]

  SEs5 <- summary(lmod5)$coefficients[,2]

  # Get intercepts and standard errors
  estIntlog <- max(intercepts5)
  estIntSE <- SEs5[which(intercepts5 == max(intercepts5))][1] # S.E. in log space

  # Convert estimate and CI to real values
  intercepts <- plogis(c(estIntlog,
                         estIntlog - (1.96 * estIntSE),
                         estIntlog + (1.96 * estIntSE))
  )
  names(intercepts) <- c("intercept.est",
                         "intercept.est.lo",
                         "intercept.est.up"
  )


  # rank intercepts (ranking is used to select SNPs colours
  # on the plot & associated lines)
  iranks <- rank(intercepts5, ties.method = "random")
  # create a data.frame associating the name of the SNP with its rank
  rankDF <- data.frame(Position=substr(names(iranks), 9, 20), rank=iranks)
  # use the merge function to add these values the appropriate row in goodDat
  goodDat <- merge(goodDat, rankDF, by="Position")

  # Put the raw data points for selected SNPs onto the rainbow plot
  with(goodDat,{
    maxx <- max(c(10, xnqlogis))
    miny <- min(c(-10, max(intercepts5)))
    plot(ylog~xnqlogis, col = rainbow(max(rank)*1.4)[rank],
         pch=1, xlim=c(0,maxx), ylim=c(miny,0),
         xlab="(Un-mapped data / mapped)",
         ylab="Allele frequency",
         main=title,
         xaxt = 'n',
         yaxt = 'n')
    }
  ) # with goodDat


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



    # Add annotation to the plot
    estDep <- plogis(-max(goodDat$xnqlogis))
    text(c(2, 7), c(0,0), c(
      paste0("Intercept est: ", signif(intercepts[1]*100, 2), "%\n",
             "(", signif(intercepts[2]*100, 2), "%-",
             signif(intercepts[3]*100, 2), "%)"),
      paste0("Mapping depth est: ", signif(estDep*100, 2), "%")
    )

    )

    # Add the trend lines to the rainbow plot
    sapply(1:length(intercepts5), function(x){
      abline(intercepts5[x], 1, col = rainbow(max(iranks)*1.4, alpha = 0.2)[iranks][x])
    })


  numgoodloci <- length(goodLoci)

  if (printout) {
    cat('Intercept based on ', numgoodloci, 'SNP loci \n')
    cat('Estimate: ', signif(intercepts[1],3), '\n')
    cat('Confindence Interval: ',
                signif(intercepts[2],3),
                '-',
                signif(intercepts[3],3),
                '\n'
                )
    cat('Mapping depth estimate: ',signif(estDep,3), '\n')
    cat('Function call \n', deparse(funcCall))
  }

  # Return the estimated values as an (initially) invisible object for further use
  invisible(list(intercepts = intercepts,
              depth.est = estDep,
              num.loci = numgoodloci,
              func.params = funcCall)
         )

}


# Get data online ####
# human
#download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/human/transformedData.csv",
#              destfile = "human.csv")
#humanDF <- read.table("human.csv", stringsAsFactors = F)

#hopper
#download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/grasshopper/transformedData.csv",
#              destfile = "hopper.csv")
#hopperDF <- read.table("hopper.csv", stringsAsFactors = F)

#parrot
#download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/parrot/transformedData.csv",
#              destfile = "parrot.csv")
#parrotDF <- read.table("parrot.csv", stringsAsFactors = F)

# Alternatively, read in data from repo ####
# setwd("~/git_repos/pseudogene_quantification/")
# # human data set
# humanDF <- read.table("data/human/transformedData.csv", stringsAsFactors = F)
# # parrot data set
# parrotDF <- read.table("data/parrot/transformedData.csv", stringsAsFactors = F)
# # grasshopper data set
# hopperDF <- read.table("data/grasshopper/transformedData.csv", stringsAsFactors = F)


# rainbowPlot(humanDF, seed = 12345, title = "Human")

# rainbowPlot(hopperDF, seed = 12345, title = "Grasshopper")

# rainbowPlot(parrotDF, seed = 12345, title = "Parrot")

