# ExoEst - Quantifying nuclear pseudogenes of non-nuclear origin from population-level data

Many eucaryotic genomes have numerous insertions of exogenous DNA (such as numts, numpts, and insertions of *Wolbachia* DNA). The quantification of these inserts can be challenging in absence of a high-quality genome assembly. ExoEst bypasses the need for such an assembly. Using population-level low-pass (genome skimming) data, it is possible to statisstically estimate the nuclear abundance of exogenous insertions. 


## Software requirements

Once the data is in the correct format, you need:
* `R` (we used version 3.6.1 (2019-07-05) -- "Action of the Toes")
* the R package `lme4` in order to construct mixed-effect models (we used version 1.1-26)

This repository includes the data sets we anaysed in the format you'll require for analysis.

If you are planning to replicate our results starting from the read data, these are the tools we used:
* `bwa` (0.7.17)
* `samtools` (1.9)
* `picard` tools (2.23.8)
* `freebayes` (1.2.0)
* GNU `parallel` (20170422)
