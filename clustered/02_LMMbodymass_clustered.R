# clustered data
# LMM size displacement testing for body mass with bioclimatic variables as covariates

# packages
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
sturnfull <- read.csv(here("data","flwbiocl.csv"),stringsAsFactors = FALSE) %>% filter(!is.na(W)) %>% filter(!is.na(sex))
occPatterns <- read.csv(here("data","spOccPatternsClust.csv"),stringsAsFactors = FALSE) %>% select(locality,sympatry)

# join to characterize localities
sturnfull <- left_join(sturnfull,occPatterns)

# units for predictors predictors
## PET to cm/month and divide temperature to original value in degrees C
sturnfull <- sturnfull %>% mutate(PETseasonality=PETseasonality/100, AnnualTmin=AnnualTmin/10)
#### refactor
sturnfull$species <- factor(sturnfull$species)
sturnfull$sympatry <- factor(sturnfull$sympatry)
sturnfull$sex <- factor(sturnfull$sex)

### model fitting 
# Full model with interaction between species and x-patry
# fit model (full)
mod1A<-lme(W~species+species:sympatry+sex+AnnualTmin+PETseasonality-1,random=~1|locality, data = sturnfull)

# variance explained
r2beta(mod1A, method = 'nsj')
# model summary
summary(mod1A)
# view effects
plot(allEffects(mod1A))

# for custom species*sympatry effect plot
BMeffectsDF <- as.data.frame(effect("species*sympatry",mod1A)) %>% mutate(measurement="bodyMass")
# write to disk
write.csv(BMeffectsDF,here("lmmEffects","clustered","Bmass.csv"),row.names = FALSE)

# anova
aovOut <- Anova(mod1A,type=2) #car package (analysis of deviance)
tidyAovOut <- aovOut %>% tidy %>%  
  rename(Chisq="statistic") %>% mutate(measurement="body mass") %>% select(measurement,everything())
# write to disk
write.csv(tidyAovOut,here("lmmEffects","clustered","aovBM.csv"),row.names = FALSE)

# fixed effects output for reporting
tidyFixed <- tidy(mod1A,effects = c("fixed")) %>% 
  select(-statistic,-p.value)
  
# confidence intervals
CIs <- intervals(mod1A)
CIs <- round(data.frame(CIs$fixed),3)

# bind estimates and CIs
bodymassLMM <- bind_cols(tidyFixed,CIs)

# calculate ICC 
varests <- as.numeric(VarCorr(mod1A)[1:2])
#vector of variance estimates
ICC <- varests[1]/sum(varests) # computes ICC

# put ICC estimate in first row of table
bodymassLMM$ICC[1] <- ICC
# reorder cols
bodymassLMM <-  bodymassLMM %>% mutate(body.measurement="body mass") %>% select(body.measurement,everything())
# write to disk
write.csv(bodymassLMM,here("lmmEffects","clustered","Woutput.csv"),row.names = FALSE)


## reduced model
# LMM with no term for interaction between species and patry
mod1A.lm1<-lme(W~species+sympatry+sex+AnnualTmin+PETseasonality-1, random=~1|locality, data = sturnfull)

# residuals from reduced model
redresW<- residuals(mod1A.lm1)

#Get model means (estimated marginal means), no interaction
W.rg1<- ref_grid(mod1A.lm1)
# results averaged over the levels of sex, but adjusted and comparable
Wsum <- summary(emmeans(W.rg1,"species","sympatry"))

# object for each predicted value
# for ease of reading
parvboth <- Wsum %>% filter(species == "parvidens" & sympatry == "both") %>% pluck("emmean")
parvallop <- Wsum %>% filter(species == "parvidens" & sympatry == "allop") %>% pluck("emmean")  
hondboth <- Wsum %>% filter(species == "hondurensis" & sympatry == "both") %>% pluck("emmean")
hondallop <- Wsum %>% filter(species == "hondurensis" & sympatry == "allop") %>% pluck("emmean")

## Adding observed residuals to no-interaction model means and get observed difference between Dsym and Dallo

# bind data and residuals
sturnAllRes<-cbind(sturnfull, redresW)
# group and summarize residual IS values
cdata <-  sturnAllRes %>% group_by(sympatry,species) %>% summarize(mean = mean(redresW, na.rm=TRUE))
# difference in trait means in sympatry
Dsym1<-((cdata[4,3]+parvboth) - (cdata[3,3] +hondboth))
Dallo1<-((cdata[2,3]+parvallop)-(cdata[1,3]+ hondallop)) 
(DDiff<- Dsym1$mean-Dallo1$mean)

# Permute residuals to whether DDiff is larger than expected by chance 
# number of permutations
nreps <- 9999
# setup vectors to put in new average residual values for sympatric and allopatric values for each species and divergence. 
parvmeansympatric<-numeric(nreps) 
parvmeanallopatric<- numeric(nreps)
hondmeansympatric<-numeric(nreps)
hondmeanallopatric<-numeric(nreps)
divergence<-numeric(nreps)

for (i in 2:nreps) {
  newW <- sample(redresW, nrow(sturnfull))
  dataset2<- cbind(sturnfull, newW)
  cdataRR <- dataset2 %>% group_by(sympatry,species) %>% summarize(mean = mean(newW, na.rm=TRUE))
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
round(c(Dsym1$mean,Dallo1$mean,DDiff,probD),3)















