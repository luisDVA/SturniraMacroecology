# NO LOCALITIES WITH N=1
# LMM size displacement testing for skulls with bioclimatic variables as covariates

# Requires
library(emmeans)
library(car)
library(dplyr)
library(purrr)
library(effects)
library(nlme)
library(r2glmm)
library(here)
library(broom)
source(here("analysis code","Baur_isosize_fn.R"))

# read the data
# working only with complete cases
batskulls <- read.csv(here("data","cranBiocl.csv"),stringsAsFactors = FALSE) %>% 
  filter(complete.cases(.)) 
occPatterns <- read.csv(here("data","spOccPatterns.csv"),stringsAsFactors = FALSE,encoding = "UTF-8") %>% select(locality,sympatry)

# remove localities with only one specimen
batskulls <-  batskulls %>% group_by(locality) %>% mutate(specimensLoc=n()) %>% filter(specimensLoc>1) %>% ungroup() %>% select(-specimensLoc)

# measurements only
Skullsmat <-  batskulls %>% select(BRH:MSF)

###################################################
#### calculate isosizes
batskulls$isosizes <- isosize(Skullsmat)[,1]

# join to characterize localities
batskulls <- left_join(batskulls,occPatterns)


#### refactor
batskulls$species <- factor(batskulls$species)
batskulls$sympatry <- factor(batskulls$sympatry)
batskulls$sex <- factor(batskulls$sex)
# modify predictors
## PET to cm/month and divide temperature to original value in degrees C
batskulls <-  batskulls%>% mutate(PETseasonality=PETseasonality/100, AnnualTmin=AnnualTmin/10)

### model fitting 
#Run full model with interaction between species and x-patry
############
# fit model (full)
mod1A<-lme(isosizes~species+species:sympatry+sex+AnnualTmin+PETseasonality-1,random=~1|locality, data = batskulls)
# for custom species*sympatry effect plot
ISeffectsDF <- as.data.frame(effect("species*sympatry",mod1A)) %>% mutate(measurement="isosize")
# write to disk
write.csv(ISeffectsDF,here("lmmEffects","filtered","IS.csv"),row.names = FALSE)
# anova
aovOut <- Anova(mod1A,type=2) #car package (analysis of deviance)
tidyAovOut <- aovOut %>% tidy %>% rename(Chisq="statistic") %>% mutate(measurement="skull isosize") %>% select(measurement,everything())
# write to disk
write.csv(tidyAovOut,here("lmmEffects","filtered","aovIS.csv"),row.names = FALSE)

# fixed effects output for reporting
tidyFixed <- tidy(mod1A,effects = c("fixed"))  %>% select(-statistic,-p.value)
# confidence intervals
CIs <- intervals(mod1A)
CIs <- data.frame(CIs$fixed)

# cbind
isosizeLMM <- bind_cols(tidyFixed,CIs)

# calculate ICC 
varests <- as.numeric(VarCorr(mod1A)[1:2])
#vector of variance estimates
ICC <- varests[1]/sum(varests) # computes ICC

# put ICC estimate in first row of table
isosizeLMM$ICC[1] <- ICC
# reorder cols
isosizeLMM <-  isosizeLMM %>% mutate(body.measurement="skull isosize") %>% select(body.measurement,everything())
# write to disk
write.csv(isosizeLMM,here("lmmEffects","filtered","ISoutput.csv"),row.names = FALSE)

#####################################
##### residual randomization 
## reduced model
# LMM with no term for interaction between species and patry
mod1A.lm1<-lme(isosizes~species+sympatry+sex+AnnualTmin+PETseasonality-1, random=~1|locality, data = batskulls)

# residuals from reduced model
redresIso<- residuals(mod1A.lm1)

#Get model means (estimated marginal means), no interaction
Iso.rg1<- ref_grid(mod1A.lm1)
# results averaged over the levels of sex, but adjusted and comparable
Isosum <- summary(emmeans(Iso.rg1,"species","sympatry"))

# object for each predicted value
# for ease of reading
parvboth <- Isosum %>% filter(species == "parvidens" & sympatry == "both") %>% pluck("emmean")
parvallop <- Isosum %>% filter(species == "parvidens" & sympatry == "allop") %>% pluck("emmean")  
hondboth <- Isosum %>% filter(species == "hondurensis" & sympatry == "both") %>% pluck("emmean")
hondallop <- Isosum %>% filter(species == "hondurensis" & sympatry == "allop") %>% pluck("emmean")

## Add observed residuals to no-interaction model means and get observed difference between Dsym and Dallo
# bind data and residuals
sturnAllRes<-cbind(batskulls, redresIso)

# group and summarize residual IS values
cdata <-  sturnAllRes %>% group_by(sympatry,species) %>% summarize(mean = mean(redresIso, na.rm=TRUE))
# difference in trait means in sympatry
Dsym1<-((cdata[4,3]+parvboth) - (cdata[3,3] +hondboth))
Dallo1<-((cdata[2,3]+parvallop)-(cdata[1,3]+ hondallop)) 
DDiff<- Dsym1$mean-Dallo1$mean

#Permute residuals to whether DDiff is larger than expected by chance 
#number of permutations
nreps <- 9999
#setup directories to put in new average residual values for sympatric and allopatric values for each species and divergence. 
parvmeansympatric<-numeric(nreps) 
parvmeanallopatric<- numeric(nreps)
hondmeansympatric<-numeric(nreps)
hondmeanallopatric<-numeric(nreps)
divergence<-numeric(nreps)

for (i in 2:nreps) {
  newIso <- sample(redresIso, nrow(batskulls))
  dataset2<- cbind(batskulls, newIso)
  cdataRR <- dataset2 %>% group_by(sympatry,species) %>% summarize(mean = mean(newIso, na.rm=TRUE))
  parvmeanallopatric[i]<-cdataRR[1,3]
  hondmeanallopatric[i]<-cdataRR[2,3]
  parvmeansympatric[i]<-cdataRR[3,3]
  hondmeansympatric[i]<-cdataRR[4,3]
  Dsym<-((cdataRR[4,3]+parvboth) - (cdataRR[3,3] +hondboth))
  Dallo<-((cdataRR[2,3]+parvallop)-(cdataRR[1,3]+hondallop)) 
  divergence[i]<- Dsym-Dallo
}

# vector of divergences from permutation
divergence <- simplify2array(divergence)

# one tailed Random permutation test 
probD <- (sum(divergence > DDiff) + 1) / (nreps+1)

# for reporting
resPermTest <- round(c(Dsym1$mean,Dallo1$mean,DDiff,probD),3) %>%(t) %>%  as_tibble
names(resPermTest) <- c("Dsym","Dallo","DDiff","pval")
resPermOut <- resPermTest %>% mutate(measurement="Skull isosize") %>% select(measurement,everything())
# write to disk
write.csv(resPermOut,here("lmmEffects","filtered","ISDDiffs.csv"),row.names = FALSE)

rm(list=ls())
