# NO LOCALITIES WITH N=1
# LMM size displacement testing for forearm length with bioclimatic variables as covariates
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

# read the data
# measurements and environment
sturnfull <- read.csv(here("data","flwbiocl.csv"),stringsAsFactors = FALSE, encoding = "UTF-8") %>% filter(!is.na(sex))
occPatterns <- read.csv(here("data","spOccPatterns.csv"),stringsAsFactors = FALSE,encoding = "UTF-8") %>% select(locality,sympatry)

# remove localities with only one specimen
sturnfull<-  sturnfull %>% group_by(locality) %>% mutate(specimensLoc=n()) %>% filter(specimensLoc>1) %>% ungroup() %>% select(-specimensLoc)

# join to characterize localities
sturnfull <- left_join(sturnfull,occPatterns)

#### refactor
sturnfull$species <- factor(sturnfull$species)
sturnfull$sympatry <- factor(sturnfull$sympatry)
sturnfull$sex <- factor(sturnfull$sex)
# modify predictors
## PET to cm/month and divide temperature to original value in degrees C
sturnfull <- sturnfull %>% mutate(PETseasonality=PETseasonality/100, AnnualTmin=AnnualTmin/10)

### model fitting 
#Run full model with interaction between species and occurrence
############

# fit model (full)
mod1A<-lme(FL~species+species:sympatry+sex+AnnualTmin+PETseasonality-1,random=~1|locality, data = sturnfull)
# for custom species*sympatry effect plot
FLeffectsDF <- as.data.frame(effect("species*sympatry",mod1A)) %>% mutate(measurement="forearm")
# write to disk
write.csv(FLeffectsDF,here("lmmEffects","filtered","FL.csv"),row.names = FALSE)
# anova
aovOut <- Anova(mod1A,type=2) #car package (analysis of deviance)
tidyAovOut <- aovOut %>% tidy %>% rename(Chisq="statistic") %>% 
  mutate(measurement="forearm length") %>% select(measurement,everything())
# write to disk
write.csv(tidyAovOut,here("lmmEffects","filtered","aovFL.csv"),row.names = FALSE)

# fixed effects output for reporting
tidyFixed <- tidy(mod1A,effects = c("fixed")) %>% select(-statistic,-p.value)
# confidence intervals
CIs <- intervals(mod1A)
CIs <- data.frame(CIs$fixed)

# cbind
forearmLMM <- bind_cols(tidyFixed,CIs)

# calculate ICC 
varests <- as.numeric(VarCorr(mod1A)[1:2])
#vector of variance estimates
ICC <- varests[1]/sum(varests) # computes ICC

# put ICC estimate in first row of table
forearmLMM$ICC[1] <- ICC
# reorder cols
forearmLMM <-  forearmLMM %>% mutate(body.measurement="forearm length") %>% select(body.measurement,everything())
# write to disk
write.csv(forearmLMM,here("lmmEffects","filtered","FLoutput.csv"),row.names = FALSE)

#####################################
##### residual randomization 
## reduced model
# LMM with no term for interaction between species and patry
mod1A.lm1<-lme(FL~species+sympatry+sex+AnnualTmin+PETseasonality-1, random=~1|locality, data = sturnfull)

# residuals from reduced model
redresFL<- residuals(mod1A.lm1)

#Get model means (estimated marginal means), no interaction
FL.rg1<- ref_grid(mod1A.lm1)
# results averaged over the levels of sex, but adjusted and comparable
FLsum <- summary(emmeans(FL.rg1,"species","sympatry"))

# object for each predicted value
# for ease of reading
parvboth <- FLsum %>% filter(species == "parvidens" & sympatry == "both") %>% pluck("emmean")
parvallop <- FLsum %>% filter(species == "parvidens" & sympatry == "allop") %>% pluck("emmean")  
hondboth <- FLsum %>% filter(species == "hondurensis" & sympatry == "both") %>% pluck("emmean")
hondallop <- FLsum %>% filter(species == "hondurensis" & sympatry == "allop") %>% pluck("emmean")

## Add observed residuals to no-interaction model means and get observed difference between Dsym and Dallo
# bind data and residuals
sturnAllRes<-cbind(sturnfull, redresFL)

# group and summarize residual IS values
cdata <-  sturnAllRes %>% group_by(sympatry,species) %>% summarize(mean = mean(redresFL, na.rm=TRUE))
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
  newFL <- sample(redresFL, nrow(sturnfull))
  dataset2<- cbind(sturnfull, newFL)
  cdataRR <- dataset2 %>% group_by(sympatry,species) %>% summarize(mean = mean(newFL, na.rm=TRUE))
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
resPermOut <- resPermTest %>% mutate(measurement="Forearm length") %>% select(measurement,everything())
# write to disk
write.csv(resPermOut,here("lmmEffects","filtered","FLDDiffs.csv"),row.names = FALSE)
rm(list=ls())
