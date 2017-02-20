### paiwise tests, cranial data, imputed datasets

#load libraries
library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)
source("modified_BESTmcmc.R")

#read the data
batSkulls <- read.csv("cranialF.csv", stringsAsFactors = F)
batSkullsMat <- batSkulls %>% dplyr::select(-cat, -sp, -locality, -collector, -sex)

## Mulitple Imputation
library(mice)

# using partial mean matching
batSkullsimp <- mice(batSkullsMat, m=10, method = "pmm")

# completed DFs
completedDFs <- list()

for (i in 1:10){ 
  completedDFs[[i]] <-  complete(batSkullsimp, i)
}

#### selecting variables for MRA analyses and setting up arguments
species <- batSkulls$sp
n <- nrow(batSkulls)


####################################################
#####  PCA in isometry free shape space
####################################################
source("Baur_MRA_functions.R")
spca <- list()
isosizes <- list()

# loop through
for (i in 1:10){
spca[[i]] <- shapePCA(completedDFs[[i]], rpc=3)
isosizes[[i]] <- as.numeric(spca[[i]]$isosize)
}
# put into DF
## V1 to V10 are each a completed dataset
isosizesDF <- as.data.frame(simplify2array(isosizes))
# bind with locality data
MRAdf <- data.frame(isosizesDF, species)
isosizeLocalities <- bind_cols(MRAdf, batSkulls) %>% select(-sp)

# find co-occurring records and allopatric records
# counting species per locality, selecting localities that have both, adding a label column and keeping unique rows

coOcurring <- 
  isosizeLocalities %>% group_by(locality)%>% 
   summarise(distinctSps=n_distinct(species)) %>% 
   filter(distinctSps>1) %>% mutate(sympatry="both") 

# merge full dataset with sympatry data
mergedWSymp <- merge(isosizeLocalities, coOcurring, by="locality", all.x=T)

# label allopatric records 
mergedWSymp$sympatry <- ifelse(is.na(mergedWSymp$sympatry), "allop", mergedWSymp$sympatry)

# La Noria and Venustiano Carranza are special cases (both species present, only hondurensis measured)
## changing manually

mergedWSymp$sympatry <- ifelse(mergedWSymp$locality=="La Noria", "both", mergedWSymp$sympatry )
mergedWSymp$sympatry <- ifelse(mergedWSymp$locality=="Venustiano Carranza", "both", mergedWSymp$sympatry)

# add column with species and occurence pattern (symp/allop)
mergedWSymp$spbiog <- paste(mergedWSymp$species, mergedWSymp$sympatry, sep = "")

# sympatric and allopatric plots

sturniraScatterplot <- 
  ggplot(mergedWSymp, aes(y=V1, x=spbiog)) +
  geom_point(aes(fill=species),pch=21,color="black",position=position_jitter(width=0.2, height=0.05)) +
  xlab("")+ylab("isosize")+
  stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, 
               geom="crossbar", width=0.3,size=0.4)+
  scale_x_discrete(labels=c(paste("S. hondurensis (A)","\n","n=",dplyr::count(mergedWSymp,spbiog)[[2]][1]),
                            paste("S. hondurensis (S)","\n","n=",dplyr::count(mergedWSymp,spbiog)[[2]][2]),
                            paste("S. parvidens (A)","\n","n=",dplyr::count(mergedWSymp,spbiog)[[2]][3]),
                            paste("S. parvidens (S)","\n","n=",dplyr::count(mergedWSymp,spbiog)[[2]][4])))+
  scale_fill_manual(values = c("#99C2CD","#FFD271"))+
  guides(fill=FALSE)


########## FOR BEST

# Split by species and co occurrence 
hondSym <- mergedWSymp %>% filter(species=="hondurensis",sympatry=="both")
hondAll <- mergedWSymp %>% filter(species=="hondurensis",sympatry=="allop")
ParvSym <- mergedWSymp %>% filter(species=="parvidens",sympatry=="both")
ParvAll <- mergedWSymp %>% filter(species=="parvidens",sympatry=="allop")

# run models Hondurensis (intraspecific)
# run once for each imputed dataset
library(BEST)
BestHondSymAllv1 <- BESTmcmc_modified(hondSym$V1, hondAll$V1)
BestHondSymAllv2 <- BESTmcmc_modified(hondSym$V2, hondAll$V2)
BestHondSymAllv3 <- BESTmcmc_modified(hondSym$V3, hondAll$V3)
BestHondSymAllv4 <- BESTmcmc_modified(hondSym$V4, hondAll$V4)
BestHondSymAllv5 <- BESTmcmc_modified(hondSym$V5, hondAll$V5)
BestHondSymAllv6 <- BESTmcmc_modified(hondSym$V6, hondAll$V6)
BestHondSymAllv7 <- BESTmcmc_modified(hondSym$V7, hondAll$V7)
BestHondSymAllv8 <- BESTmcmc_modified(hondSym$V8, hondAll$V8)
BestHondSymAllv9 <- BESTmcmc_modified(hondSym$V9, hondAll$V9)
BestHondSymAllv10 <- BESTmcmc_modified(hondSym$V10, hondAll$V10)

# pool the MCMC chains
library(runjags)
combinedHondSymAll <- combine.mcmc(list(BestHondSymAllv1$samples, BestHondSymAllv2$samples,
                      BestHondSymAllv3$samples, BestHondSymAllv4$samples,
                      BestHondSymAllv5$samples, BestHondSymAllv6$samples,
                      BestHondSymAllv7$samples, BestHondSymAllv8$samples,
                      BestHondSymAllv9$samples, BestHondSymAllv10$samples))

# wrangle the pooled chains into a results DF
mcmcChain = as.matrix(combinedHondSymAll)
  colnames(mcmcChain) <- c("mu1", "mu2", "nu", "sigma1", 
                           "sigma2")
mcmcDF <- as.data.frame(mcmcChain)
class(mcmcDF) <- c("BEST", class(mcmcDF))

# summarize for reporting
summary(mcmcDF,credmass=0.99)
hondSymAllimpCombined <- mcmcDF
# write to disk (optional)
saveRDS(object = hondSymAllimpCombined, file="BESThondSymAllimpComb.rds")

### repeat for parvidens (intraspecific)
BestParvSymAllv1 <- BESTmcmc_modified(ParvSym$V1, ParvAll$V1)
BestParvSymAllv2 <- BESTmcmc_modified(ParvSym$V2, ParvAll$V2)
BestParvSymAllv3 <- BESTmcmc_modified(ParvSym$V3, ParvAll$V3)
BestParvSymAllv4 <- BESTmcmc_modified(ParvSym$V4, ParvAll$V4)
BestParvSymAllv5 <- BESTmcmc_modified(ParvSym$V5, ParvAll$V5)
BestParvSymAllv6 <- BESTmcmc_modified(ParvSym$V6, ParvAll$V6)
BestParvSymAllv7 <- BESTmcmc_modified(ParvSym$V7, ParvAll$V7)
BestParvSymAllv8 <- BESTmcmc_modified(ParvSym$V8, ParvAll$V8)
BestParvSymAllv9 <- BESTmcmc_modified(ParvSym$V9, ParvAll$V9)
BestParvSymAllv10 <- BESTmcmc_modified(ParvSym$V10, ParvAll$V10)

# pool chains
combinedParvSymAll <- combine.mcmc(list(BestParvSymAllv1$samples, BestParvSymAllv2$samples,
                                        BestParvSymAllv3$samples, BestParvSymAllv4$samples,
                                        BestParvSymAllv5$samples, BestParvSymAllv6$samples,
                                        BestParvSymAllv7$samples, BestParvSymAllv8$samples,
                                        BestParvSymAllv9$samples, BestParvSymAllv10$samples))

# wrangle the pooled chain into a results DF
mcmcChain = as.matrix(combinedParvSymAll)
colnames(mcmcChain) <- c("mu1", "mu2", "nu", "sigma1", 
                         "sigma2")
mcmcDF <- as.data.frame(mcmcChain)
class(mcmcDF) <- c("BEST", class(mcmcDF))
summary(mcmcDF,credMass=0.99)
round(summary(mcmcDF,credMass=0.99),3)
ParvSymAllimpCombined <- mcmcDF

# write to disk (optional)
saveRDS(object = ParvSymAllimpCombined, file="BESTParvSymAllimpComb.rds")

#overall differences (interspecific)
hondfullImp <- MRAdf %>% filter(species=="hondurensis")
parvfullImp <- MRAdf %>% filter(species=="parvidens")

BestOverallSpImp1 <- BESTmcmc_modified(hondfullImp$V1, parvfullImp$V1)
BestOverallSpImp2 <- BESTmcmc_modified(hondfullImp$V2, parvfullImp$V2)
BestOverallSpImp3 <- BESTmcmc_modified(hondfullImp$V3, parvfullImp$V3)
BestOverallSpImp4 <- BESTmcmc_modified(hondfullImp$V4, parvfullImp$V4)
BestOverallSpImp5 <- BESTmcmc_modified(hondfullImp$V5, parvfullImp$V5)
BestOverallSpImp6 <- BESTmcmc_modified(hondfullImp$V6, parvfullImp$V6)
BestOverallSpImp7 <- BESTmcmc_modified(hondfullImp$V7, parvfullImp$V7)
BestOverallSpImp8 <- BESTmcmc_modified(hondfullImp$V8, parvfullImp$V8)
BestOverallSpImp9 <- BESTmcmc_modified(hondfullImp$V9, parvfullImp$V9)
BestOverallSpImp10 <- BESTmcmc_modified(hondfullImp$V10, parvfullImp$V10)

# pool chains
combinedOverallSp <- combine.mcmc(list(BestOverallSpImp1$samples, BestOverallSpImp2$samples,
                                       BestOverallSpImp3$samples, BestOverallSpImp4$samples,
                                       BestOverallSpImp5$samples, BestOverallSpImp6$samples,
                                       BestOverallSpImp7$samples, BestOverallSpImp8$samples,
                                       BestOverallSpImp9$samples, BestOverallSpImp10$samples))

# wrangle the pooled chains into a results DF
mcmcChain = as.matrix(combinedOverallSp)
colnames(mcmcChain) <- c("mu1", "mu2", "nu", "sigma1", 
                         "sigma2")
mcmcDF <- as.data.frame(mcmcChain)
class(mcmcDF) <- c("BEST", class(mcmcDF))
round(summary(mcmcDF,credMass=0.99),3)
combinedOverallSp <- mcmcDF

# write to disk
saveRDS(object = combinedOverallSp, file="combinedOverallSp.rds")

### further interspecific tests
# both species Sympatric/Sympatric
BestHondParvSymv1 <- BESTmcmc_modified(hondSym$V1, ParvSym$V1)
BestHondParvSymv2 <- BESTmcmc_modified(hondSym$V2, ParvSym$V2)
BestHondParvSymv3 <- BESTmcmc_modified(hondSym$V3, ParvSym$V3)
BestHondParvSymv4 <- BESTmcmc_modified(hondSym$V4, ParvSym$V4)
BestHondParvSymv5 <- BESTmcmc_modified(hondSym$V5, ParvSym$V5)
BestHondParvSymv6 <- BESTmcmc_modified(hondSym$V6, ParvSym$V6)
BestHondParvSymv7 <- BESTmcmc_modified(hondSym$V7, ParvSym$V7)
BestHondParvSymv8 <- BESTmcmc_modified(hondSym$V8, ParvSym$V8)
BestHondParvSymv9 <- BESTmcmc_modified(hondSym$V9, ParvSym$V9)
BestHondParvSymv10 <- BESTmcmc_modified(hondSym$V10, ParvSym$V10)

# pool chains
combinedHondParvSym <- combine.mcmc(list(BestHondParvSymv1$samples,BestHondParvSymv2$samples,
                                        BestHondParvSymv3$samples,BestHondParvSymv4$samples,
                                        BestHondParvSymv5$samples,BestHondParvSymv6$samples,
                                        BestHondParvSymv7$samples,BestHondParvSymv8$samples,
                                        BestHondParvSymv9$samples,BestHondParvSymv10$samples))
# repeat earlier step for getting summary
mcmcChain = as.matrix(combinedHondParvSym)
colnames(mcmcChain) <- c("mu1", "mu2", "nu", "sigma1", 
                         "sigma2")
mcmcDF <- as.data.frame(mcmcChain)
class(mcmcDF) <- c("BEST", class(mcmcDF))
round(summary(mcmcDF,credMass=0.99),3)
HondParvSymimpCombined <- mcmcDF

# write to disk
saveRDS(object = HondParvSymimpCombined, file="HondParvSymimpCombined.rds")

################
# both species Allopatric/Allopatric
BestHondParvAllv1 <- BESTmcmc_modified(ParvAll$V1, hondAll$V1)
BestHondParvAllv2 <- BESTmcmc_modified(ParvAll$V2, hondAll$V2)
BestHondParvAllv3 <- BESTmcmc_modified(ParvAll$V3, hondAll$V3)
BestHondParvAllv4 <- BESTmcmc_modified(ParvAll$V4, hondAll$V4)
BestHondParvAllv5 <- BESTmcmc_modified(ParvAll$V5, hondAll$V5)
BestHondParvAllv6 <- BESTmcmc_modified(ParvAll$V6, hondAll$V6)
BestHondParvAllv7 <- BESTmcmc_modified(ParvAll$V7, hondAll$V7)
BestHondParvAllv8 <- BESTmcmc_modified(ParvAll$V8, hondAll$V8)
BestHondParvAllv9 <- BESTmcmc_modified(ParvAll$V9, hondAll$V9)
BestHondParvAllv10 <- BESTmcmc_modified(ParvAll$V10, hondAll$V10)

# pool chains
combinedHondParvAll<- combine.mcmc(list(BestHondParvAllv1$samples, BestHondParvAllv2$samples,
                                         BestHondParvAllv3$samples, BestHondParvAllv4$samples,
                                         BestHondParvAllv5$samples, BestHondParvAllv6$samples,
                                         BestHondParvAllv7$samples, BestHondParvAllv8$samples,
                                         BestHondParvAllv9$samples, BestHondParvAllv10$samples))
# combine and summarize
mcmcChain = as.matrix(combinedHondParvAll)
colnames(mcmcChain) <- c("mu1", "mu2", "nu", "sigma1", 
                         "sigma2")
mcmcDF <- as.data.frame(mcmcChain)
class(mcmcDF) <- c("BEST", class(mcmcDF))
round(summary(mcmcDF,credMass=0.99),3)
HondParvAllimpCombined <- mcmcDF
# write to disk
saveRDS(object = HondParvAllimpCombined, file="HondParvAllimpCombined.rds")
