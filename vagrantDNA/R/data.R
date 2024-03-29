#' NUMT datasets
#'
#' These data sets give allele frequencies generated from whole-genome
#' sequencing data mapped against mitochondrial genome references.
#'
#' @format Each is a data.frame.
#' \strong{\code{[species]DF}} just contain variant calling data.
#' The columns are:
#' \describe{
#'   \item{Sample}{A character vector giving a unique name for each sample}
#'   \item{Position}{A numeric vector giving the positio of each SNP}
#'   \item{AltProp}{A numeric vector giving the proportion (of reads mapping to the exogenous genome) that carry the non-standard alleles (thought to be in the vagrant copies)}
#'   \item{DP}{A numeric vector giving the mapping depth at each site}
#'   \item{nMapped}{A numeric vector giving the number of base pairs in this sample's sequencing data that were successfully aligned to the extranuclear reference}
#'   \item{nTot}{A numeric vector giving the total number of base pairs of the sample's mapping data (ideally after quality control, filtering, read trimming, ect.)}
#'   \item{ylog}{A numeric vector giving the log of `AltProp`}
#' }
#' \strong{\code{[species]FX}} were generated from species with diverged populations. These
#' contain information on sites with fixed differences. The columns are:
#' \describe{
#'   \item{pos}{A vector giving the site IDs}
#'   \item{sample}{A vector, giving a unique name for each individual genotyped}
#'   \item{g1}{A numeric vector giving the allele count of allele 1}
#'   \item{g2}{A numeric vector giving the allele count of allele 2}
#'   \item{g3}{A numeric vector giving the allele count of allele 3}
#'   \item{g4}{A numeric vector giving the allele count of allele 4}
#'   \item{N}{A numeric vector giving the number of bp in the read data that did not map to the vargrant DNA reference in this individual}
#'   \item{M}{A numeric vector giving the number of bp in the read data that did map to the vargrant DNA reference in this individual}
#'   \item{pop}{A factor (or structure that can be coerced to a factor), of "A" and "B" denoting which population the individual belongs to}
#'   \item{A}{A numeric vector giving the number of the major allele at this site in population A}
#'   \item{B}{A numeric vector giving the number of the major allele at this site in population B}
#'   }
#'
#' @examples
#' # Generate the estimates reported in the paper
#' \dontrun{
#' ## Parrot
#' rPar <- rainbowPlot(parrotDF, seed=12345)
#' # remove multiple outliers
#' toRemove <- interceptPositionPlot(rPar, highlightOutliers=T)
#' parrotDF2 <- parrotDF[!parrotDF$Position %in% toRemove,]
#' rainbowPlot(parrotDF2, seed=12345, title = "Parrot")
#'
#' ## Human
#' rHum <- rainbowPlot(humanDF, seed=12345)
#' # automatic outlier detection would remove SNPs that look OK
#' interceptPositionPlot(rHum, highlightOutliers = T)
#' # select single outlier "by hand"
#' which.max(coef(summary(rHum$lmer.model))[,1])
#' humanDF2 <- humanDF[humanDF$Position != 310,] # by hand, better
#' rainbowPlot(humanDF2, seed=12345, title = "Human")
#'
#' ## Grasshopper
#' download.file("https://tinyurl.com/4mtrbkzc", destfile = "hopper.csv")
#' hopperDF <- read.table("hopper.csv")
#' rainbowPlot(hopperDF, seed = 12345, title = "Grasshopper")
#' }
#' humanDF
#' parrotDF
#' hopperFX
#' parrotFX
"humanDF"

#' @rdname humanDF
#' @format NULL
"parrotDF"

#' @rdname humanDF
#' @format NULL
"hopperFX"

#' @rdname humanDF
#' @format NULL
"parrotFX"
