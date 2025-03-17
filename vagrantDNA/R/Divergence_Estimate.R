# A function to generate the divergence estimate of vagrant DNA inserts.

#' Estimate vagrant insert proportion from diverged population data
#'
#' @param dat A data.frame with the following eleven columns in this order.
#' \describe{
#'   \item{pos}{A factor (or structure that can be coerced to a factor),
#'   giving the site IDs}
#'   \item{sample}{A factor (or structure that can be coerced to a factor),
#'   giving a unique name for each individual genotyped}
#'   \item{g1}{A numeric vector giving the allele count of allele 1.}
#'   \item{g2}{A numeric vector giving the allele count of allele 2.}
#'   \item{g3}{A numeric vector giving the allele count of allele 3.}
#'   \item{g4}{A numeric vector giving the allele count of allele 4.}
#'   \item{N}{A numeric vector giving the number of bp in the read data that did not map
#'    to the vargrant DNA reference in this individual.}
#'   \item{M}{A numeric vector giving the number of bp in the read data that did map
#'   to the vargrant DNA reference in this individual.}
#'   \item{pop}{A factor (or structure that can be coerced to a factor), of "A" and "B"
#'   denoting which population the individual belongs to.}
#'   \item{A}{A numeric vector giving the number of the major allele at this site
#'   in population A.}
#'   \item{B}{A numeric vector giving the number of the major allele at this site
#'   in population B.}
#'   }
#'
#' @return An (invisible) list of Ascores and Bscores corresponding to the
#' estimates per locus using A or B as the focal population.
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @export
#' @examples
#' ## Access one of the package's example data-sets (parrotFX or hopperFX)
#' data(hopperFX)

#'
#' ## plot and printout (by default) the results of running divEst on the grasshopper
#' divEst(hopperFX)
#' ##
#' ## save the estimates from each locus and population and locus
#' res1 <- divEst(hopperFX)
#' ## Inspect the individual scores
#' str(res1)
#' \dontrun{
#'
#'
#' }
divEst <- function(dat){

  siteNames <- unique(dat$pos)
  nSites <- length(siteNames)
  Ascores <- Bscores <- rep(0,nSites)

  dat$A1 <- dat$A == 1
  dat$A2 <- dat$A == 2
  dat$A3 <- dat$A == 3
  dat$A4 <- dat$A == 4
  dat$B1 <- dat$B == 1
  dat$B2 <- dat$B == 2
  dat$B3 <- dat$B == 3
  dat$B4 <- dat$B == 4
  dat$allDep <- rowSums(dat[,3:6])
  dat$aMaj <- rowSums(dat[,3:6] * dat[,12:15])
  dat$aAlt <- rowSums(dat[,3:6] * !dat[,12:15])
  dat$bMaj <- rowSums(dat[,3:6] * dat[,16:19])
  dat$bAlt <- rowSums(dat[,3:6] * !dat[,16:19])

  for (i in 1:nSites){# create matrices with counts for the A populations and B populations
    gmatA <- as.matrix(dat[dat$pos==siteNames[i] & dat$pop == 'A',3:6])
    gmatB <- as.matrix(dat[dat$pos==siteNames[i] & dat$pop == 'B',3:6])

    # create matrices identifying the non mito alleles in each population
    # and the mito allele in the opposite population (where it will be non-mito)
    # For the A population:
    notMitoAlleleA <- !as.matrix(dat[dat$pos==siteNames[i] & dat$pop == 'A',12:15])
    AmitoAlleleInB <-  as.matrix(dat[dat$pos==siteNames[i] & dat$pop == 'B',12:15])
    # For the B populations
    notMitoAlleleB <- !as.matrix(dat[dat$pos==siteNames[i] & dat$pop == 'B',16:19])
    BmitoAlleleInA <-  as.matrix(dat[dat$pos==siteNames[i] & dat$pop == 'A',16:19])
    mA <- dat[dat$pos==siteNames[i] & dat$pop == 'A',]$M
    mB <- dat[dat$pos==siteNames[i] & dat$pop == 'B',]$M
    nA <- dat[dat$pos==siteNames[i] & dat$pop == 'A',]$N
    nB <- dat[dat$pos==siteNames[i] & dat$pop == 'B',]$N

    Ascores[i] <- sum(rowSums(gmatA * notMitoAlleleA) / rowSums(gmatA) * mA) / sum(nA) +
      sum(rowSums(gmatB * AmitoAlleleInB) / rowSums(gmatB) * mB) / sum(nB)
    Bscores[i] <- sum(rowSums(gmatB * notMitoAlleleB) / rowSums(gmatB) * mB) / sum(nB) +
      sum(rowSums(gmatA * BmitoAlleleInA) / rowSums(gmatA) * mA) / sum(nA)
  }


  plot(Ascores, Bscores)
  abline(0,1)
  abline(coef(lm(Bscores ~ Ascores)), col = "Blue")
  abline(v = mean(Ascores), col='Orange')
  abline(h = mean(Bscores), col='Orange')

  print(noquote(paste("Mean of A Scores: ",
                      signif(mean(Ascores), 3),
                      " (SE ",
                      signif(sqrt(var(Ascores)/nSites),2),
                      ")"
  )))

  print(noquote(paste("Mean of B Scores: ",
                      signif(mean(Bscores), 3),
                      " (SE ",
                      signif(sqrt(var(Bscores)/nSites),2),
                      ")"
  )))
  invisible(list(
    Ascores=Ascores,
    Bscores=Bscores
  ))
}
