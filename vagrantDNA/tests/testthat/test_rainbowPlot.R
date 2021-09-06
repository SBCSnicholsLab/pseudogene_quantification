library("checkmate")

# creating a test data.frame
Sample <- rep(paste('sample', 4^(1:6), sep = ''), each = 4) # Sample names: six samples scored at four loci
pE <- rep(1/4^(1:6), each = 4) # Proportion extranuclear reads in each sample
pN <- 1-pE # Proportion of nuclear reads in each sample
Position <- rep(paste('locus', 1:4, sep = ''), 6) # Locus names: four loci recorded for six samples
pA <- rep(c(0.9, 0.8, 0.7, 0.6), 6) # frequency of alternative allele at the nuMT loci

pV <- 0.0001 # proportion of each genome that is vagrant

AltProp <- pN*pV*pA / (pN*pV + pE) # proportion of reads mapping to the extra-nuclear genome that have alternate allele
rNE <- pN*(1-pV) / (pE + pN*pV)  # ratio of reads mapping to extranuclear genome (includes vagrants) / those that do not
xnqlogis <- log(rNE)  # log and add error so regression has residual var

isd <- 1e-03 # sd per individual
rsd <- 1e-07 # residual sd
# add this variation to log(AltProp)
ylog <- log(AltProp) + rep(rnorm(6, sd = isd), each = 4) + rnorm(6*4, sd = rsd)


testDF <- data.frame(
  Sample = Sample,
  xnqlogis =xnqlogis,
  Position = Position,
  ylog = ylog,
  AltProp = AltProp
  )

# testDF gives a good fit with maximum intercept 9.00009e-05 and SE ~ 0
testlist <- rainbowPlot(testDF, nloci = 4, minSamples = 6)
test_that("The function returns appropriate objects", {
  expect_list(testlist, types = c('vector','numeric','integerish','lmerMod','call'))
})

test_that("The parameter estimates are correct for a perfect fit", {
  expect_equal(as.numeric(testlist$intercepts), rep(9.00009e-05, 3), tolerance = 0.01)
  expect_equal(testlist$depth.est,1-plogis(max(testDF$xnqlogis)))
  expect_equal(testlist$num.loci, 4)
})
