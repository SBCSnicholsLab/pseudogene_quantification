
#' Intercept-position-plot for loci seect by rainbowPlot()
#'
#' @param ff A `rainbowPlot()` fit.
#' @param main Plot caption, default is "Intercept-position-plot"
#' @param ... Arguments to be passed to `plot()`
#'
#' @export
#'
#' @examples
#'  \dontrun{
#'  humanFit <- rainbowPlot(humanDF)
#'  interceptPositionPlot(humanFit)
#'  }
interceptPositionPlot <- function(ff, main="Intercept-position-plot", ...){
  coefs.a <- as.data.frame(coef(summary(ff$lmer.model)))
  coefs.a$site.numeric <- as.numeric(substr(rownames(coefs.a), 9, 30))
  plot(Estimate ~ site.numeric,
       data=coefs.a,
       main=main,
       ...)
  l01 <- loess(Estimate ~ site.numeric,
               data=coefs.a)
  lines(predict(l01, coefs.a$site.numeric)~coefs.a$site.numeric)
}


