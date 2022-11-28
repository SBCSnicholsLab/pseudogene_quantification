



#' Boxplot of alternate (vagrant) allele counts per sample
#'
#' @param dd A dataframe, such as used to run `rainbowPlot()`
#' @param ff (optional) a `rainbowPlot()` fit
#' @param ... Arguments passed to `boxplot()`
#'
#' @export
#'
#' @examples
#' #'  \dontrun{
#'  humanFit <- rainbowPlot(humanDF)
#'  selAllPlot(humanDF) # overall alt allele counts
#'  selAllPlot(humanDF, humanFit) # alt allele count at loci selected for fit
#'  }
selAllPlot <- function(dd,ff,ylim=c(0, 10),...){
  if(hasArg(ff)){
    posString <- rownames(coef(summary(ff$lmer.model)))
    sites <- substr(posString, 9, 20)
    pDat <- dd[dd$Position %in% sites,]
    boxplot(pDat$DP * exp(pDat$ylog) ~ pDat$Sample,
            main="Alt. allele count at selected loci",
            ylim=ylim,
            ...)
    goodDatMeds <- round(as.numeric(tapply(pDat$DP * exp(pDat$ylog), pDat$Sample, function(x) boxplot.stats(x)$stats[3])))
    if(sum(goodDatMeds < 2) > length(goodDatMeds)/2) warning("The majority of individuals have a median alternate allele depth of 1. This suggests that vagrant DNA proportion or the sequencing depth may be insufficient for an accurate estimate. Consider using more WGS data if possible.")
  } else {
    boxplot(dd$DP * exp(dd$ylog) ~ dd$Sample,
            main="Alt. allele count at all loci",
            ylim=ylim,
            ...)
    cat("Plot shows alt allels counts for all loci. Supply rainbowPlot object to show only selected loci.\n")
    goodDatMeds <- round(as.numeric(tapply(dd$DP * exp(dd$ylog), dd$Sample, function(x) boxplot.stats(x)$stats[3])))
    if(sum(goodDatMeds < 2) > length(goodDatMeds)/2) warning("The majority of individuals have a median alternate allele depth of 1. This suggests that vagrant DNA proportion or the sequencing depth may be insufficient for an accurate estimate. Consider using more WGS data if possible.")
  }
}

