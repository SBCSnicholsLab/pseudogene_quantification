## code to prepare humanDF, hopperDF and ParrotDF

# human
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/human/transformedData.csv",
              destfile = "human.csv")
humanDF <- read.table("human.csv")

#hopper
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/grasshopper/transformedData.csv",
              destfile = "hopper.csv")
hopperDF <- read.table("hopper.csv")

#parrot
download.file("https://raw.githubusercontent.com/SBCSnicholsLab/pseudogene_quantification/main/data/parrot/transformedData.csv",
              destfile = "parrot.csv")
parrotDF <- read.table("parrot.csv")

usethis::use_data(humanDF, hopperDF, parrotDF, overwrite = TRUE)
