
#' Intercept-position-plot for loci selected by \code{rainbowPlot()}
#'
#' @param ff A \code{rainbowPlot()} fit.
#' @param main Plot caption, default is "Intercept-position-plot"
#' @param highlightOutliers Logical. Whether or not to highlight sites with
#' outlier intercepts according to \code{boxplot.stats}. Default \code{F}
#' @param ... Arguments to be passed to \code{plot()}
#' @return Returns invisibly a vector of positions where the intercept estimate
#' were outliers according to \code{boxplot.stats}
#'
#' @export
#'
#' @examples
#'  \dontrun{
#'  rPar <- rainbowPlot(parrotDF)
#'  toRemove <- interceptPositionPlot(rPar, highlightOutliers=T)
#'  parrotDF2 <- parrotDF[!parrotDF$Position %in% toRemove,]
#'  rainbowPlot(parrotDF2)
#'  }
interceptPositionPlot <- function(ff,
                                  highlightOutliers=F,
                                  main="Intercept-position-plot",
                                  ...){
  coefs.a <- as.data.frame(coef(summary(ff$lmer.model)))
  coefs.a$site.numeric <- as.numeric(substr(rownames(coefs.a), 9, 30))
  bounds <- boxplot.stats(coefs.a$Estimate)$stat
  coefs.a$out <- coefs.a$Estimate < bounds[1] | coefs.a$Estimate > bounds[5]
  if(highlightOutliers) {
  plot(Estimate ~ site.numeric,
       data=coefs.a,
       col=out+1,
       main=main,
       ...)
  } else {
    plot(Estimate ~ site.numeric,
         data=coefs.a,
         main=main,
         ...)
  }
  l01 <- loess(Estimate ~ site.numeric,
               data=coefs.a)
  lines(predict(l01, coefs.a$site.numeric)~coefs.a$site.numeric)
  return(invisible(coefs.a$site.numeric[coefs.a$out]))
}


