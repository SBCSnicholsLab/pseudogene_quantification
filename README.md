# vagrantDNA - Quantifying nuclear inserts of etxra-nuclear origin from population-level data

Many eukaryotic genomes contain numerous inserts of extra nucler DNA (such as numts, numpts, and insertions of *Wolbachia* DNA). The quantification of these inserts can be challenging in absence of a high-quality genome assembly. ExtraIns bypasses the need for an assembly. Using population-level low-pass (genome skimming) data, it is possible to statistically estimate the nuclear abundance of such inserts. 


## Setup and software requirements
### Setup
```
# Install devtools from CRAN
install.packages("devtools")

# Install vagrantDNA from GitHub
devtools::install_github("https://github.com/SBCSnicholsLab/pseudogene_quantification/",
                         subdir = "vagrantDNA")
# Load package
library(vagrantDNA)
```

### Generate a rainbowplot
```
download.file("t.ly/6hXO", destfile = "hopper.csv")
hopperDF <- read.table("~/hopper.csv")
res1 <- rainbowPlot(parrotDF, seed = 12345, printout = FALSE, title = "Grasshopper")
## print just the stored estimates (the first two elements of the list)
print(res1[1:2])
## Inspect the residuals of the lmer model
plot(res1$lmer.model)

```


## Reanalysis

If you are planning to replicate our results starting from the read data, these are the tools we used:
* `bwa` (0.7.17)
* `samtools` (1.9)
* `picard` tools (2.23.8)
* `freebayes` (1.2.0)
* GNU `parallel` (20170422)

The sequencing data are available from the SRA.
