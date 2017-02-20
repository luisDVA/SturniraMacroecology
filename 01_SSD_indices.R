# Sexual size dimorphism indices 
library(dplyr)

# read FL data and remove specimens with missing sex value
sturnFull <- read.csv("finalFLdataset.csv", stringsAsFactors = FALSE)
sturnFull <- sturnFull %>% filter(sex!="NA")

# split by species then by sex
hondurensis <- sturnFull %>% filter(species=="hondurensis")
hondurensisM <- hondurensis %>% filter(sex=="M")
hondurensisF <- hondurensis %>% filter(sex=="F")
parvidens <- sturnFull %>% filter(species=="parvidens")
parvidensM <- parvidens %>% filter(sex=="M")
parvidensF <- parvidens %>% filter(sex=="F")

# get SSD ratios for FL (Stephens and Wiens 2009)
# S. hondurensis
mean(hondurensisM$FL)
mean(hondurensisF$FL)
# males are larger so they go in the numerator and the SSD index is multiplied
### by minus 1 
#### rounded for reporting   
round(-1*(mean(hondurensisM$FL)/mean(hondurensisF$FL)-1),2)

# repeat for parvidens
mean(parvidensM$FL)
mean(parvidensF$FL)
# males are also larger so they go in the numerator and the SSD index is multiplied
### by minus 1 
#### rounded for reporting   
round(-1 *(mean(parvidensM$FL)/mean(parvidensF$FL)- 1),2)

