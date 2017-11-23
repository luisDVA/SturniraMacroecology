# LMM size displacement testing for skulls with bioclimatic variables as covariates
#Requires
library(emmeans)
library(car)
library(dplyr)
library(purrr)
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
#### calculate isosizes
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


#### refactor
mergedWSymp$species <- factor(mergedWSymp$species)
mergedWSymp$sympatry <- factor(mergedWSymp$sympatry)
mergedWSymp$sex <- factor(mergedWSymp$sex)

############
# fit model (full)
mod1A<-lmer(isosizes~species*sympatry+sex+AnnualTmin+PETseasonality+(1|locality), data = mergedWSymp)
summary(mod1A)
plot(allEffects(mod1A))



#####################################
###  bioclimatic covariates
##### residual randomization


### model fitting 
#Run full model with interaction between species and patry
mod1A<-lmer(isosizes~species*sympatry+AnnualTmin+PETseasonality+(1|locality), data = mergedWSymp)
# anova
Anova(mod1A, contrasts=list(elevation=contr.sum, species=contr.sum, sympatry=contr.sum), type=3)
# for reporting
round(Anova(mod1A, contrasts=list(elevation=contr.sum, species=contr.sum, sympatry=contr.sum), type=3),2)


#Create LMM with no term for interaction between species and patry
mod1A.lm1<-lmer(isosizes~species+sympatry+AnnualTmin+PETseasonality+(1|locality), data = mergedWSymp)
# anova
Anova(mod1A.lm1, contrasts=list(elevation=contr.sum, species=contr.sum, sympatry=contr.sum), type=3)

# residuals from reduced model
redresIso<- residuals(mod1A.lm1)

#Get model means (estimated marginal means), no interaction
Iso.rg1<- ref_grid(mod1A.lm1)
#summarize
Isosum<-summary(Iso.rg1)
Isosum
# object for each predicted value
# for ease of reading
parvboth <- Isosum %>% filter(species == "parvidens" & sympatry == "both") %>% pluck("prediction")
parvallop <- Isosum %>% filter(species == "parvidens" & sympatry == "allop") %>% pluck("prediction")  
hondboth <- Isosum %>% filter(species == "hondurensis" & sympatry == "both") %>% pluck("prediction")
hondallop <- Isosum %>% filter(species == "hondurensis" & sympatry == "allop") %>% pluck("prediction")

## Add observed residuals to no-interaction model means and get observed difference between Dsym and Dallo
# bind data and residuals
sturnAllRes<-cbind(mergedWSymp, redresIso)

# group and summarize residual IS values
cdata <-  sturnAllRes %>% group_by(sympatry,species) %>% summarize(mean = mean(redresIso, na.rm=TRUE))
# difference in trait means in sympatry
Dsym1<-((cdata[4,3]+parvboth) - (cdata[3,3] +hondboth))
Dallo1<-((cdata[2,3]+parvallop)-(cdata[1,3]+ hondallop)) 
DDiff<- Dsym1-Dallo1
round(Dsym1,3)
round(Dallo1,3)
round(DDiff,3)

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
  newIso <- sample(redresIso, 196)
  dataset2<- cbind(mergedWSymp, newIso)
  cdataRR <- dataset2 %>% group_by(sympatry,species) %>% summarize(mean = mean(newIso, na.rm=TRUE))
  parvmeanallopatric[i]<-cdataRR[1,3]
  hondmeanallopatric[i]<-cdataRR[2,3]
  parvmeansympatric[i]<-cdataRR[3,3]
  hondmeansympatric[i]<-cdataRR[4,3]
  Dsym<-((cdataRR[4,3]+parvboth) - (cdataRR[3,3] +hondboth))
  Dallo<-((cdataRR[2,3]+parvallop)-(cdataRR[1,3]+hondallop)) 
  divergence[i]<- Dsym-Dallo
}


#Calculate the probability that observed diverences is greater than random
probD<-length(divergence[divergence>= DDiff$mean])/nreps[1]
probD

# for reporting
round(c(Dsym1$mean,Dallo1$mean,DDiff$mean,probD),2)

#Write a file with diveregence values for each permutation in case anyone thinks you made it up
#write.csv(divergence, file="ISdivergence.csv")















