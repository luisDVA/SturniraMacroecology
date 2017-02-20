#########################
# Cranial data, isosize comparisons with complete cases

# load libraries
library(ggplot2)
library(cowplot)
library(gridExtra)
source("Baur_MRA_functions.R")
library(dplyr)

#read the data
batSkulls <- read.csv("cranialF.csv", stringsAsFactors = F) %>% filter(complete.cases(.))
batSkullsMat <- batSkulls %>% dplyr::select(-cat, -sp, -locality, -collector, -sex)

#### selecting variables for analyses and setting up arguments
species <- batSkulls$sp
n <- nrow(batSkulls)


####################################################
#####  PCA in isometry free shape space
####################################################

# run shape PCA and calculate isometric size axis
spca <- shapePCA(batSkullsMat, rpc=3)
isosizes<- as.numeric(spca$isosize)

# merge with locality data
MRAdf <- data.frame(isosizes, species)
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

######
# sympatric and allopatric plots

library(ggforce)

sturniraScatterplot <- 
ggplot(mergedWSymp, aes(y=isosizes, x=spbiog)) +
  geom_sina(aes(fill=species),pch=21,color="black",maxwidth=0.35) +
  xlab("")+ylab("isosize")+
  stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, 
               geom="crossbar", width=0.4,size=0.4)+
  scale_x_discrete(labels=c(paste("S. hondurensis (A)","\n","n=",dplyr::count(mergedWSymp,spbiog)[[2]][1]),
                            paste("S. hondurensis (S)","\n","n=",dplyr::count(mergedWSymp,spbiog)[[2]][2]),
                            paste("S. parvidens (A)","\n","n=",dplyr::count(mergedWSymp,spbiog)[[2]][3]),
                            paste("S. parvidens (S)","\n","n=",dplyr::count(mergedWSymp,spbiog)[[2]][4])))+
  scale_fill_manual(values = c("#335EAD","#EEBB33"))+
  guides(fill=FALSE)


########## FOR BEST

# Split by species and co occurrence 
hondSym <- mergedWSymp %>% filter(species=="hondurensis", sympatry=="both")
hondAll <- mergedWSymp %>% filter(species=="hondurensis", sympatry=="allop")
ParvSym <- mergedWSymp %>% filter(species=="parvidens", sympatry=="both")
ParvAll <- mergedWSymp %>% filter(species=="parvidens", sympatry=="allop")

library(BEST)

# run models S. hondurensis (intraspecific)
BestHondSymAllcomp <- BESTmcmc(hondSym$isosizes, hondAll$isosizes)
# for reporting results
round(summary(BestHondSymAllcomp, credmass=0.99),3)
# write to disk
saveRDS(object = BestHondSymAllcomp, file="BestHondSymAllcomp.rds")

### repeat for parvidens (intraspecific)
BestParvSymAllcomp <- BESTmcmc(ParvSym$isosizes,ParvAll$isosizes)
BestParvSymAllcomp <- readRDS("BesParvS")
# for reporting results
round(summary(BestParvSymAllcomp,credMass=0.99),3)
# write to disk
saveRDS(object = BestParvSymAllcomp, file="BestParvSymAllcomp.rds")

# for figure
# save summary objects
summHondBEST <- summary(BestHondSymAllcomp,credMass = 0.99) %>% as.data.frame()
summParvBEST <- summary(BestParvSymAllcomp,credMass = 0.99) %>% as.data.frame()
# save plot objects
HondDiffs <- plot(BestHondSymAllcomp,credMass = 0.99)
ParvDiffs <- plot(BestParvSymAllcomp,credMass = 0.99)

# make plot objects into DFs
hondDiffsDF <- data.frame(counts=HondDiffs$counts,
                          dens=HondDiffs$density,
                          mids=HondDiffs$mids)

parvDiffsDF <- data.frame(counts=ParvDiffs$counts,
                          dens=ParvDiffs$density,
                          mids=ParvDiffs$mids)

# plot objects
hondDiffsgrob <- 
  ggplot(hondDiffsDF,aes(x=mids,y=dens))+geom_bar(stat="identity",fill="#335EAD")+
  geom_segment(aes(x = summHondBEST["muDiff","HDIlo"],
                   y = 0, xend = summHondBEST["muDiff","HDIup"], 
                   yend = 0),size=2.5)+
  xlab(expression(mu[1] - mu[2]))+ylab("density")+
  geom_segment(aes(x=0,y=0,xend=0,yend=max(hondDiffsDF$dens)*.75),linetype=2)+
  annotate("text",x=0,y=max(hondDiffsDF$dens)*.79,
           label="0")

parvDiffsgrob <- 
  ggplot(parvDiffsDF,aes(x=mids,y=dens))+geom_bar(stat="identity",fill="#EEBB33")+
  geom_segment(aes(x = summParvBEST["muDiff","HDIlo"],
                   y = 0, xend = summParvBEST["muDiff","HDIup"], 
                   yend = 0),size=2.5)+
  xlab(expression(mu[1] - mu[2]))+ylab("density")+
  geom_segment(aes(x=0,y=0,xend=0,yend=max(parvDiffsDF$dens)*.88),linetype=2)+
  annotate("text",x=0,y=max(parvDiffsDF$dens)*.91,
           label="0")

# put the plots together
grid.arrange(sturniraScatterplot, arrangeGrob(hondDiffsgrob, parvDiffsgrob, ncol=2), ncol=1, heights=c(1,0.52))

#overall differences (interspecific)
hondfull <- mergedWSymp %>% filter(species=="hondurensis")
Parvfull <- mergedWSymp %>% filter(species=="parvidens")
# interspecific comparison
BestOverallSp <- BESTmcmc(hondfull$isosizes, Parvfull$isosizes)
# optional, write object to disk
saveRDS(object = BestOverallSp, file="bestOverall.rds")
# for reporting
round(summary(BestOverallSp, credMass = 0.99),3)


# more models Symp/Symp (interspecific)
BestHondParvSymcomp <- BESTmcmc(hondSym$isosizes, ParvSym$isosizes)
# for reporting results
round(summary(BestHondParvSymcomp, credMass=0.99),3)

# more models A/A (interspecific)
BestHondParvAllcomp <- BESTmcmc(hondAll$isosizes, ParvAll$isosizes)
# for reporting results
round(summary(BestHondParvAllcomp, credMass=0.99),3)
