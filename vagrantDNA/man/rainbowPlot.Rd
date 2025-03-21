% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rainbow_Plot.R
\name{rainbowPlot}
\alias{rainbowPlot}
\title{Rainbow Plot.}
\usage{
rainbowPlot(
  data,
  nloci = 400,
  minWt = 0.01,
  maxFreq = 0.7,
  minSamples = 10,
  seed,
  title = "",
  printout = TRUE,
  correctForDepth = F,
  weigh = T,
  extraNucLen = 16000
)
}
\arguments{
\item{data}{A data.frame with at least the following columns.
\describe{
  \item{Sample}{A factor (or structure that can be coerced to a factor),
  giving a unique name for each sample.}
  \item{Position}{A factor (or structure that can be coerced to a factor),
  giving a unique name for each SNP location.}
  \item{AltProp}{A numeric vector giving the proportion (of reads mapping to
  the exogenous genome) that carry the non-standard alleles (thought to be in
  the vagrant copies)}
  \item{DP}{A numeric vector giving the mapping depth at each site}
  \item{nMapped}{A numeric vector giving the number of base pairs in this
  sample's sequencing data that were successfully aligned to the extranuclear
  reference}
  \item{nTot}{A numeric vector giving the total number of base pairs of the sample's
  mapping data (ideally after quality control, filtering, read trimming, ect.)}
  \item{ylog}{A numeric vector giving the log(AltProp)}
  }}

\item{nloci}{The number of loci to be selected for the analysis. Default, 400.}

\item{minWt}{The minimum average allele frequency of SNPs to be included in the analysis.
Default, 0.01.}

\item{maxFreq}{The maximum allele frequency of an individual observation to be included
in the analysis. Default, 0.7.}

\item{minSamples}{The minimum number of samples in which a SNP should be called in order to
be included in the analysis.}

\item{seed}{Random number seed.}

\item{title}{User-supplied title for the rainbow plot.}

\item{printout}{If printout is TRUE, the function prints the estimates. Default, `TRUE`.}

\item{correctForDepth}{Logical, whether of not to correct for uneven insertion rates along the extranuclear sequence. Not usually required. Default is `FALSE`.}

\item{weigh}{Logical. Whether or not to select loci with high allele frequencies. Default is `TRUE`.}

\item{extraNucLen}{Numeric. Length of the extranuclear genome reference. Only required when `correctForDepth` is set to `TRUE`. Default is 16000.}
}
\value{
An invisible list with the following elements.
\describe{
 \item{$intercepts}{A vector giving the intercept estimate,
 the lower and upper 95\% confidence interval bounds.}
 \item{$depth.est}{A crude upper limit; the minimum over samples
 of the proportion of reads mapping to the exogenous genome.
 Some part (or all) of the CI would be expected to be below this value.}
 \item{$num.loci}{The number of loci remaining after filtering,
 which were used to obtain the intercept estimates.}
 \item{$lmer.model}{An lmer object which contains the details of the mixed effects model fitted to identify
 the largest intercept and its standard error.
 In the case of unusual or problematic results it may be useful to inspect the residuals, see examples.}
 \item{$func.params}{A record of the function call}
}
}
\description{
Generate a plot to estimate the proportion of the nuclear genome
made up of a particular vagrant sequence. The function produces the plots
intercept estimate, and mapping depth estimate described by
Becher & Nichols (2023) Molecular Ecology Resources.
}
\examples{
## Access one of the package's example data-sets (parrotDF or humanDF)
data(parrotDF)
##
## n.b. hopperDF is too large to be included in the package's data
## but can be accessed as follows
 \dontrun{
 download.file("https://tinyurl.com/4mtrbkzc", destfile = "hopper.csv")
 hopperDF <- read.table("hopper.csv")}

## (the t.ly/ link is to the CSV file at
## https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/
## main/data/grasshopper/transformedData.csv)


## plot and printout (by default) the results of running rainbowPlot on the parrot data.
rainbowPlot(parrotDF, seed = 12345, title = "Parrot")
##
## plot without printing the results and store results in parrotFit
parrotFit <- rainbowPlot(parrotDF, seed = 12345, printout = FALSE, title = "Parrot")
## Inspect relationship between intercepts and SNP location
interceptPositionPlot(parrotFit)
## Inspect mapping depth of nuclear insert allele across samples
selAllPlot(parrotDF, parrotFit)

}
