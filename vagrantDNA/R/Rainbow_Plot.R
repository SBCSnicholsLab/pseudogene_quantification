################################################################################
# A method to estimate the nuclear frequency of exogenous pseudogenes from     #
#   population-level low-pass sequencing data                                  #
################################################################################



#' Rainbow Plot.
#'
#' Generate a plot to estimate the proportion of the nuclear genome
#' made up of a particular vagrant sequence. The function produces the plots
#' intercept estimate, and mapping depth estimate described by
#' Becher & Nichols (2021) citation.
#'
#' @importFrom lme4 lmer
#' @importFrom grDevices boxplot.stats rainbow
#' @importFrom graphics abline axis text
#' @importFrom methods hasArg
#' @importFrom stats coef complete.cases lm plogis var
#' @importFrom utils install.packages
#'
#' @param data A data.frame with at least the following columns.
#' \describe{
#'   \item{Sample}{A factor (or structure that can be coerced to a factor),
#'   giving a unique name for each sample.}
#'   \item{Position}{A factor (or structure that can be coerced to a factor),
#'   giving a unique name for each SNP location.}
#'   \item{AltProp}{A numeric vector giving the proportion (of reads mapping to
#'   the exogenous genome) that carry the non-standard alleles (thought to be in
#'   the vagrant copies)}
#'   \item{DP}{A numeric vector giving the mapping depth at each site}
#'   \item{nMapped}{A numeric vector giving the number of base pairs in this
#'   sample's sequencing data that were successfully aligned to the extranuclear
#'   reference}
#'   \item{nTot}{A numeric vector giving the total number of base pairs of the sample's
#'   mapping data (ideally after quality control, filtering, read trimming, ect.)}
#'   \item{ylog}{A numeric vector giving the log(AltProp)}
#'   }
#' @param nloci The number of loci to be selected for the analysis. Default, 400.
#' @param minWt The minimum average allele frequency of SNPs to be included in the analysis.
#' Default, 0.01.
#' @param maxFreq The maximum allele frequency of an individual observation to be included
#' in the analysis. Default, 0.7.
#' @param minSamples The minimum number of samples in which a SNP should be called in order to
#' be included in the analysis.
#' @param seed Random number seed.
#' @param title User-supplied title for the rainbow plot.
#' @param printout If printout is TRUE, the function prints the estimates. Default, `TRUE`.
#' @param correctForDepth Logical, whether of not to correct for uneven insertion rates along the extranuclear sequence. Not usually required. Default is `FALSE`.
#' @param weigh Logical. Whether or not to select loci with high allele frequencies. Default is `TRUE`.
#' @param extraNucLen Numeric. Length of the extranuclear genome reference. Only required when `correctForDepth` is set to `TRUE`. Default is 16000.
#' @return An invisible list with the following elements.
#' \describe{
#'  \item{$intercepts}{A vector giving the intercept estimate,
#'  the lower and upper 95\% confidence interval bounds.}
#'  \item{$depth.est}{A crude upper limit; the minimum over samples
#'  of the proportion of reads mapping to the exogenous genome.
#'  Some part (or all) of the CI would be expected to be below this value.}
#'  \item{$num.loci}{The number of loci remaining after filtering,
#'  which were used to obtain the intercept estimates.}
#'  \item{$lmer.model}{An lmer object which contains the details of the mixed effects model fitted to identify
#'  the largest intercept and its standard error.
#'  In the case of unusual or problematic results it may be useful to inspect the residuals, see examples.}
#'  \item{$func.params}{A record of the function call}
#' }
#'
#' @examples
#' ## Access one of the package's example data-sets (parrotDF or humanDF)
#' data(parrotDF)
#' ##
#' ## n.b. hopperDF is too large to be included in the package's data
#' ## but can be accessed as follows
#'  \dontrun{
#'  download.file("https://tinyurl.com/4mtrbkzc", destfile = "hopper.csv")
#'  hopperDF <- read.table("hopper.csv")}
#'
#' ## (the t.ly/ link is to the cvs file at
#' ## https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/
#' ## main/data/grasshopper/transformedData.csv)
#'
#'
#' ## plot and printout (by default) the results of running rainbowPlot on the parrot data.
#' rainbowPlot(parrotDF, seed = 12345, title = "Parrot")
#' ##
#' ## plot without printing the results and store results in parrotFit
#' parrotFit <- rainbowPlot(parrotDF, seed = 12345, printout = FALSE, title = "Parrot")
#' ## Inspect relationship between intercepts and SNP location
#' interceptPositionPlot(parrotFit)
#' ## Inspect mapping depth of nuclear insert allele across samples
#' selAllPlot(parrotDF, parrotFit)
#'
#' @export
rainbowPlot <- function(data,
                      nloci = 400,
                      minWt = 0.01,
                      maxFreq = 0.7,
                      minSamples = 10,
                      seed,
                      title = "",
                      printout = TRUE,
                      correctForDepth=F,
                      weigh=T,
                      extraNucLen=16000
                      ){
  if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")
  # force Position and Sample to become factors
  data$Position <- as.factor(data$Position)
  data$Sample <- as.factor(data$Sample)
  if(correctForDepth){
    cat("Option correctForDepth is set to true. Make sure to specify the length or the extranulcear genome reference using the argument extraNucLen!\n")
    p <- exp(data$ylog)
    nma <- p * data$DP
    ma <- (1-p) * data$DP
    rd <- data$DP
    ard <- data$nMapped/extraNucLen
    data$ylog <- log(nma*ard/rd/((ard/rd)*nma+ma))
  }

  # check the data.frame has appropriate columns
  if( any( c( !is.vector(data$AltProp),
              !is.factor(data$Position),
              !is.vector(data$ylog),
              !is.factor(data$Sample),
              !is.vector(data$nMapped),
              !is.vector(data$nTot),
              !is.vector(data$DP)
              )
           )
      ) stop("Supply a dataframe containing a factors Position & Sample (or vectors which can be coerced to a factor), \n
              \t plus vectors AltProp, ylog & xnqlogisBP")

  # Greate required columns
  data$mappingrate <- data$nMapped / data$nTot
  data$xnqlogis <- -qlogis(data$mappingrate)
  data <- data[!is.na(data$xnqlogis),]
  # Save the function call
  funcCall <- sys.call()

  # Remove the rows of data with NAs
  goodDat <- data[complete.cases(data),]

  # Remove abnormally high values
  goodDat <- goodDat[goodDat$AltProp < maxFreq,]


  # Find how many individuals each SNP has been scored in
  tt <- table(goodDat$Position)
  # get the names of SNPs occurring in more than minSamples
  goodNames <- names(tt)[tt>=minSamples]

  # Exclude loci found in too few samples
  goodDat <- goodDat[goodDat$Position %in% goodNames,]

  # set random number seed if requested (ie. if seed is specified in the function call)
  if(hasArg(seed)) set.seed(seed)

  # Find the average allele frequency at each Position
  wweights <- tapply(goodDat$AltProp, factor(goodDat$Position), mean)
  # exclude loci with average frequency less than minWt
  wweights[wweights < minWt] <- 0

  if (sum(wweights>0) < nloci) stop("Your filtering steps have resulted in too few loci for sampling consider \n
                                    \t decreasing nloci or minSamples,
                                    \t or perhaps decreasing maxFreq or minWt")

  # sample a set of high frequency loci (sampling proportional to the freq)
  if(weigh){
    loci <- sample(names(wweights), nloci, prob = wweights)
  } else {
    loci <- sample(names(wweights), nloci)
  }


  # Calculate the raw slopes for each locus (to identify any outlying values)
  slopes <- rep(0, nloci)
  for (i in 1:nloci) {
    df <- goodDat[goodDat$Position == loci[i],]
    #if(perBpDep) {
    #slopes[i] <- coef( lm(ylog ~ xnqlogisBP, data = df) )[2]
    #} else {
      slopes[i] <- coef( lm(ylog ~ xnqlogis, data = df) )[2]
    #}
    }

  # Round slopes to avoid problems with rounding errors
  slopes <- round(slopes, 5)

  # Use the boxplot.stats function to identifying SNPs with outlying slopes
  # Hard filtering excludes upper and lower quartiles
  # otherwise the boxplot default criterion is used (familiar as the end of whiskers)
  bpsSlopes <- boxplot.stats(slopes)$stats
  #if (filterHard) {
   # rogueSNPs <- (slopes > bpsSlopes[4] | slopes < bpsSlopes[2])}
  #else {
    rogueSNPs <- (slopes > bpsSlopes[5] | slopes < bpsSlopes[1])
    #}
  goodLoci <- loci[!rogueSNPs]

  # remove the rogue loci from the data.frame
  goodDat <- goodDat[goodDat$Position %in% goodLoci,]

  # get median from boxplot.stats (discards outliers)
  goodDatMeds <- round(as.numeric(tapply(goodDat$DP * exp(goodDat$ylog), goodDat$Sample, function(x) boxplot.stats(x)$stats[3])))

  if(sum(goodDatMeds < 2) > length(goodDatMeds)/2) warning("The majority of individuals have a median alternate allele depth of 1. This suggests that vagrant DNA proportion or the sequencing depth may be insufficient for an accurate estimate. Consider using more WGS data if possible.")
  # fit a 1:1 line plus intercept plus sample random effect
  #if(perBpDep) {
  #lmod5 <- lmer(ylog ~ 0 + Position + (1 | Sample),
  #              offset = xnqlogisBP,
  #              data = goodDat)
  #} else {
    lmod5 <- lmer(ylog ~ 0 + Position + (1 | Sample),
                  offset = xnqlogis,
                  data = goodDat)
  #}


  intercepts5 <- summary(lmod5)$coefficients[,1]

  SEs5 <- summary(lmod5)$coefficients[,2]

  # distances between intercepts
  # intDiffs <- diff(sort(intercepts5, decreasing = T))
  # thresh <- median(intDiffs) * 5
  # lastAccepted <- min(which(intDiffs > thresh))


  # Get intercepts and standard errors
  estIntlog <- max(intercepts5)
  #estIntlog <- sort(intercepts5, decreasing = T)[lastAccepted]
  #estIntSE <- SEs5[which(intercepts5 == sort(intercepts5, decreasing = T)[lastAccepted])][1] # S.E. in log space
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
  # if(perBpDep){
  #   maxx <- max( c(10, goodDat$xnqlogisBP) )
  #   miny <- min( c(-10, max(intercepts5)) )
  #   plot(goodDat$ylog ~ goodDat$xnqlogisBP,
  #        col = rainbow(max(goodDat$rank)*1.4)[goodDat$rank],
  #        pch=1, xlim=c(0,maxx), ylim=c(miny,0),
  #        xlab="(Un-mapped data / mapped)",
  #        ylab="Allele frequency",
  #        main=title,
  #        xaxt = 'n',
  #        yaxt = 'n'
  #        )
  # } else {
    maxx <- max( c(10, goodDat$xnqlogis) )
    miny <- min( c(-10, max(intercepts5)) )
    plot(goodDat$ylog ~ goodDat$xnqlogis,
         col = rainbow(max(goodDat$rank)*1.4)[goodDat$rank],
         pch=1, xlim=c(0,maxx), ylim=c(miny,0),
         xlab="(Un-mapped data / mapped)",
         ylab="Allele frequency",
         main=title,
         xaxt = 'n',
         yaxt = 'n'
    )
#}

  # add the axes

    axis(1, at = log(2^(0:7*2)), labels = 2^(0:7*2))
    abline(v = log(2^(0:7*2)), col = "grey", lty=3)
    ll <- expression("1", "10"^-1, "10"^-2, "10"^-3, "10"^-4, "10"^-5, "10"^-6)
    axis(2, at = log(10^(0:-6)), labels = ll)
    abline(h = log(10^(0:-6)), col = "grey", lty=3)

    # Add key lines
    abline(h=max(intercepts5))
    #abline(h=sort(intercepts5, decreasing = T)[lastAccepted])
    # if(perBpDep){
    # abline(v=max(goodDat$xnqlogisBP))
    # } else {
      abline(v=max(goodDat$xnqlogis))
    #}
    abline(max(intercepts5), 1,  lty=2)
    #abline(sort(intercepts5, decreasing = T)[lastAccepted], 1,  lty=2)



    # Add annotation to the plot
    #if(perBpDep){
    #estDep <- plogis(-max(goodDat$xnqlogisBP))
    #} else {
      estDep <- plogis(-max(goodDat$xnqlogis))
    #}
    text(c(2, 2), c(0,-1), c(
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
    cat('Confidence Interval: ',
                signif(intercepts[2],3),
                '-',
                signif(intercepts[3],3),
                '\n'
                )
    cat('Mapping depth estimate: ',signif(estDep,3), '\n')
    cat('Function call \n', deparse(funcCall), '\n')
  }

  if(estDep < intercepts[1]) warning("The intercept estimate is higher than the mapping-depth estimate. This fit is likely to be unreliable. Consider using more sequencing data per individual if possible.")

  # Return the estimated values as an (initially) invisible object for further use
  invisible(list(intercepts = intercepts,
              depth.est = estDep,
              num.loci = numgoodloci,
              lmer.model = lmod5,
              func.params = funcCall
              )
         )

}

