\name{hopperDF}
\alias{hopperDF}
\docType{data}
\title{The grasshopper SNP dataframe}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("hopperDF")}
\format{
  A data frame with 401360 observations on the following 7 variables.
  \describe{
    \item{\code{Sample}}{a character vector (that can be coerced to a factor), giving a unique name for each sample.}
    \item{\code{Position}}{a character vector (that can be coerced to a factor), giving a unique name for each sample.}
    \item{\code{AltProp}}{A numberic vector giving the proportion (of reads mapping to the exogenous genome) that carry the non-standard alleles thought to be in the vagrant copies}
    \item{\code{mappingrate}}{a numeric vector giving the proportion of reads that map to the exogenous genome (intermediate value in the calculation of xnqlogis not used by the rainbowPlot function)}
    \item{\code{ylog}}{a numeric vector, log(AltProp)}
    \item{\code{xnqlogis}}{a numeric vector, A numeric vector giving log(m/N); where m is the number of reads mapping to the exogenous genome and N is the remaining reads.}
    \item{\code{xlog}}{a numeric vector, a vector used in intermediate calculations, not used by the rainbowplot function}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(hopperDF)
## maybe str(hopperDF) ; plot(hopperDF) ...
}
\keyword{datasets}
