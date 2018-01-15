# LMM size displacement testing for skulls with bioclimatic variables as covariates
#Requires
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
batskulls <- read.csv(here("data","cranBiocl.csv"),stringsAsFactors = FALSE, fileEncoding = "UTF-8")
# occurrence in the localities
occPatterns <- read.csv(here("data","spOccPatterns.csv"),stringsAsFactors = FALSE,encoding = "UTF-8") %>% select(locality,sympatry)

# join to characterize localities
batskulls <- left_join(batskulls,occPatterns)

#### refactor
batskulls$species <- factor(batskulls$species)
batskulls$sympatry <- factor(batskulls$sympatry)
batskulls$sex <- factor(batskulls$sex)

# units for predictors predictors
## PET to cm/month and divide temperature to original value in degrees C
batskulls <- batskulls %>% mutate(PETseasonality=PETseasonality/100, AnnualTmin=AnnualTmin/10)

### model fitting 
# Full model with interaction between species and x-patry
mod1A<-lme(CBL~species*sympatry+sex+AnnualTmin+PETseasonality,random=~1|locality, data = batskulls)
# model r2
r2beta(mod1A, method = 'nsj')
# view effects
plot(allEffects(mod1A))
# for custom species*sympatry effect plot
headEffectsDF <- as.data.frame(effect("species*sympatry",mod1A)) %>% mutate(measurement="CBL")
# write to disk
write.csv(headEffectsDF,here("lmmEffects","CBL.csv"),row.names = FALSE)
# anova
aovOut <- Anova(mod1A,type=2) #car package (Wald Chisq)
tidyAovOut <- aovOut %>% tidy %>% 
  rename(Chisq="statistic") %>% mutate(measurement="head length") %>% select(measurement,everything())
# write to disk
write.csv(tidyAovOut,here("lmmEffects","aovCBL.csv"),row.names = FALSE)

# fixed effects output for reporting
tidyFixed <- tidy(mod1A,effects = c("fixed")) %>% select(-statistic,-p.value)
# confidence intervals
CIs <- intervals(mod1A)
CIs <- data.frame(CIs$fixed)

# bind estimates and CIs
headlengthLMM <- bind_cols(tidyFixed,CIs)

# calculate ICC 
varests <- as.numeric(VarCorr(mod1A)[1:2])
#vector of variance estimates
ICC <- varests[1]/sum(varests) # computes ICC

# put ICC estimate in first row of table
headlengthLMM$ICC[1] <- ICC
# reorder cols
headlengthLMM <-  headlengthLMM %>% mutate(body.measurement="head length") %>% select(body.measurement,everything())
# write to disk
write.csv(headlengthLMM,here("lmmEffects","CBLoutput.csv"),row.names = FALSE)


## reduced model
# LMM with no term for interaction between species and occurrence
mod1A.lm1<-lme(CBL~species+sympatry+sex+AnnualTmin+PETseasonality, random=~1|locality, data = batskulls)

# residuals from reduced model
redresCBL<- residuals(mod1A.lm1)

# Get model means (estimated marginal means), no interaction
CBL.rg1<- ref_grid(mod1A.lm1)
# results averaged over the levels of sex, but adjusted and comparable
CBLsum <- summary(emmeans(CBL.rg1,"species","sympatry"))

# object for each predicted value
# for ease of reading
parvboth <- CBLsum %>% filter(species == "parvidens" & sympatry == "both") %>% pluck("emmean")
parvallop <- CBLsum %>% filter(species == "parvidens" & sympatry == "allop") %>% pluck("emmean")  
hondboth <- CBLsum %>% filter(species == "hondurensis" & sympatry == "both") %>% pluck("emmean")
hondallop <- CBLsum %>% filter(species == "hondurensis" & sympatry == "allop") %>% pluck("emmean")

## Adding observed residuals to no-interaction model means and get observed difference between Dsym and Dallo

# bind data and residuals
sturnAllRes<-cbind(batskulls, redresCBL)

# group and summarize residual IS values
cdata <-  sturnAllRes %>% group_by(sympatry,species) %>% summarize(mean = mean(redresCBL, na.rm=TRUE))
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
  newCBL <- sample(redresCBL, nrow(batskulls))
  dataset2<- cbind(batskulls, newCBL)
  cdataRR <- dataset2 %>% group_by(sympatry,species) %>% summarize(mean = mean(newCBL, na.rm=TRUE))
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

# proportion of random divergences greater than observed
probD<-length(divergence[divergence>= DDiff])/nreps[1]

# for reporting
round(c(Dsym1$mean,Dallo1$mean,DDiff,probD),3)

# view
hist(divergence)
abline(v=-DDiff)














