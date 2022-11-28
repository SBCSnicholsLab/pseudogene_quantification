## code to prepare humanDF, hopperDF and ParrotDF

# human
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/human02/transformedData.csv",
              destfile = "human.csv")
humanDF <- read.table("human.csv")


#hopper
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/grasshopper/transformedData.csv",
              destfile = "hopper.csv")
hopperDF <- read.table("hopper.csv")

download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/grasshopper/hopperFixed.csv",
              destfile = "hopperFixed.csv")
hopperFX <- read.table("hopperFixed.csv")


#parrot
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/parrot/transformedData.csv",
              destfile = "parrot.csv")
parrotDF <- read.table("parrot.csv")

download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/parrot/parrotFixed.csv",
              destfile = "parrotFixed.csv")
parrotFX <- read.table("parrotFixed.csv")


usethis::use_data(humanDF, parrotDF, parrotFX, hopperFX, overwrite = TRUE, compress = "bzip2")
