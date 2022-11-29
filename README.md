# vagrantDNA - Quantifying nuclear inserts of etxra-nuclear origin from population-level data

Many eukaryotic genomes contain numerous inserts of extra nucler DNA (such as nuMts [mitochondrial inserts], numPs [plasmid], and insertions of *Wolbachia* DNA). We use the term 'vagrant DNA' to cover all these different types of inserted sequence. The assembly of vagrant DNA is problematic, so it is difficult to estimate how much is in the genome. 

Instead the package vagrantDNA estimates the proportion of the nuclear genome is made up of a particular vagrant sequence without the need for genome assembly. The the method is designed to use cheap low-pass (genome skimming) data. 


## Setup and software requirements
### Setup
```
# Install devtools from CRAN
install.packages("devtools")

# Install vagrantDNA from GitHub
devtools::install_url("https://github.com/SBCSnicholsLab/pseudogene_quantification/releases/download/v1.1.0/vagrantDNA_1.1.0.tar.gz")
# Load package
library(vagrantDNA)
```

### Generate a rainbowplot
```
download.file("https://tinyurl.com/4mtrbkzc", destfile = "hopper.csv")
hopperDF <- read.table("hopper.csv")
hopperFit <- rainbowPlot(hopperDF, seed = 12345, printout = FALSE, title = "Grasshopper")
```


## Reanalysis

If you are planning to replicate our results starting from the read data, these are the tools we used:
* `bwa` (0.7.17)
* `samtools` (1.9)
* `picard` tools (2.23.8)
* `freebayes` (1.2.0)
* GNU `parallel` (20170422)

The sequencing data are available from the SRA.
