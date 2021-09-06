\name{parrotDF}
\alias{parrotDF}
\docType{data}
\title{The parrot SNP dataframe}
\description{A dataset required designed for the function rainbowPlot.
Each row refers to a particular SNP locus and a particular sample.
The column Sample gives a unique identifier for each sample, and Position for each SNP locus.
For the reads that map to the extranuclear genome in question (mitochondria in this case),
each row gives the proportion of reads at that SNP which carry the alleles thought to be in the vagrant sequences,
i.e. the alleles which are not found in that sample's mitochondria .  The raw frequency of the putative vagrant alleles is AltProp, and its log is ylog. Note some the vagrant sequences,
inserted into the parrot genome may carry the extranuclear allele: that is what makes estimating the frequency of vagrant
sequences problematic. Each row also gives the log(N/m) where m is count of reads that map to the mitochondrion, and N is the count of those that do not.  The latter will usually be from the nuclear genome, but could include contamination.}
\usage{data("parrotDF")}
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
\details{details of the pipeline to produce this dataframe can be found at https://github.com/SBCSnicholsLab/pseudogene_quantification/tree/main/data/parrot}
\source{
from Becher and Nichols 2022
}
\references{background information in }

\examples{
data(parrotDF); str(parrotDF); plot(parrotDF)
}
\keyword{datasets}
