\name{parrotFX}
\alias{parrotFX}
\docType{data}
\title{
Grasshopper Divergence Estimate Dataset
}
\description{
A dataset designed for the function divEst, which calculated the proportion of the nuclear genome comprising vagrant DNA. In the case of this dataset divEst will calculation the proportion of the nuclear genome of the parrot Psephotellus varius which is NUMTs, making use of the fact that individuals from different populations(populations A and B) carry different mitochondrial genomes, and the alleles at the SNPs listed here differ between populations A and B.
}
\usage{data("parrotFX")}
\format{
  A data frame with 5772 observations on the following 11 variables.
  \describe{
    \item{\code{pos}}{A factor (or structure that can be coerced to a factor), giving the site IDs.}
    \item{\code{sample}}{A factor (or structure that can be coerced to a factor), giving a unique name for each individual genotyped.}
    \item{\code{g1}}{A numeric vector giving the allele count of allele 1.}
    \item{\code{g2}}{A numeric vector giving the allele count of allele 2.}
    \item{\code{g3}}{A numeric vector giving the allele count of allele 3.}
    \item{\code{g4}}{A numeric vector giving the allel count of allele 4.}
    \item{\code{N}}{A numeric vector giving the number of bp in the read data that did NOT map to the mtDNA reference in this individual.}
    \item{\code{M}}{A numeric vector giving the number of bp in the read data that did map to the mtDNA reference in this individual.}
    \item{\code{pop}}{A factor (or structure that can be coerced to a factor), of "A" and "B" denoting which population the individual belongs to.}
    \item{\code{A}}{A numeric vector giving the number of the major allele at this site in population A.}
    \item{\code{B}}{A numeric vector giving the number of the major allele at this site in population B.}
  }
}
\details{
details of the pipeline to produce this dataframe can be found at https://github.com/SBCSnicholsLab/pseudogene_quantification/tree/main/data/parrot
}
\source{
https://github.com/SBCSnicholsLab/pseudogene_quantification/blob/main/data/grasshopper/hopperFixed.csv
}
\references{
from Becher H and Nichols RA (2022). Assembly-free quantification of vagrant DNA inserts
}
\examples{
data(parrotFX)
maybe str(parrotFX)
divEst(parrotFX)
}
\keyword{datasets}
