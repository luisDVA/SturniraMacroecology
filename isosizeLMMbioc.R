# LMM size displacement testing for skulls with bioclimatic variables as covariates
#Requires
#library(lsmeans)
#library(car)
library(dplyr)
library(purrr)
library(ggplot2)
library(effects)
library(lme4)
source("Baur_isosize_fn.R")

# read the data
# working only with complete cases
batSkulls <- read.csv("cranBiocl.csv",stringsAsFactors = FALSE) %>% filter(complete.cases(.))

# split by species
# parvidens
parvall <- batSkulls %>% filter(sp=="parvidens")
parvSkullsmat <-  batSkulls %>% filter(sp=="parvidens") %>% select(BRH:MSF)
# hondurensis
hondall <- batSkulls %>% filter(sp=="hondurensis")
hondSkullsmat <-  batSkulls %>% filter(sp=="hondurensis") %>% select(BRH:MSF)

###################################################
#### isosizes

parvall$isosizes <- isosize(parvSkullsmat)[,1]
hondall$isosizes <- isosize(hondSkullsmat)[,1]

# bind into new DF
batSkullsF <- bind_rows(parvall,hondall)
batSkullsF <-  batSkullsF %>% rename(species=sp)

# find co-occurring records and allopatric records
# counting species per locality, selecting localities that have both, adding a label column and keeping unique rows

coOcurring <- 
  batSkullsF %>% group_by(locality)%>% 
  summarise(distinctSps=n_distinct(species)) %>% 
  filter(distinctSps>1) %>% dplyr::mutate(sympatry="both") 

# merge full dataset with sympatry data
mergedWSymp <- merge(batSkullsF,coOcurring,by="locality",all.x=T)

# label allopatric records 
mergedWSymp$sympatry <- ifelse(is.na(mergedWSymp$sympatry),"allop",mergedWSymp$sympatry)

# La Noria is a special case (both species present, only hondurensis measured)
## changing manually
mergedWSymp$sympatry <- ifelse(mergedWSymp$locality=="La Noria","both",mergedWSymp$sympatry )
# add column with species and occurence pattern (symp/allop)
mergedWSymp$spbiog <- paste(mergedWSymp$species,mergedWSymp$sympatry,sep = "")

#
#### refactor
mergedWSymp$species <- factor(mergedWSymp$species)
mergedWSymp$sympatry <- factor(mergedWSymp$sympatry)
mergedWSymp$sex <- factor(mergedWSymp$sex)

############
# fit model
mod1A<-lmer(isosizes~species*sympatry+sex+b1+b12+(1|locality), data = mergedWSymp)
summary(mod1A)
plot(allEffects(mod1A))