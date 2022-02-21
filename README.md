# vagrantDNA - Quantifying nuclear inserts of etxra-nuclear origin from population-level data

Many eukaryotic genomes contain numerous inserts of extra nucler DNA (such as nuMts [mitochondrial inserts], numPs [plasmid], and insertions of *Wolbachia* DNA). We use the term 'vagrant DNA' to cover all these different types of inserted sequence. The assembly of vagrant DNA is problematic, so it is difficult to estimate how much is in the genome. 

Instead the package vagrantDNA estimates the proportion of the nuclear genome is made up of a particular vagrant sequence without the need for genome assembly. The the method is designed to use cheap low-pass (genome skimming) data. 


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
