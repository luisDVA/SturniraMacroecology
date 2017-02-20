# pairwise tests, Forearm Length

# load packages (install first if needed)

library(dplyr)
library(ggplot2)
library(cowplot)
library(BEST)
library(gridExtra)
library(ggforce)

# read dataset
sturnFull <- read.csv("finalFLdataset.csv",stringsAsFactors = FALSE)
  
# find co-occurring records and allopatric records
# counting species per locality, selecting localities that have both, adding a label column and keeping unique rows

coOcurring <- 
  sturnFull %>% group_by(locality)%>% 
  summarise(distinctSps=n_distinct(species)) %>% 
  filter(distinctSps>1) %>% mutate(sympatry="both") 


# merge full dataset with sympatry data
mergedWSymp <- merge(sturnFull, coOcurring, by="locality", all.x=T)

# label allopatric records 
mergedWSymp$sympatry <- ifelse(is.na(mergedWSymp$sympatry), "allop", mergedWSymp$sympatry)

# add column with species and occurence pattern (symp/allop)
mergedWSymp$spbiog <- paste(mergedWSymp$species, mergedWSymp$sympatry, sep = "")


# sympatric and allopatric plots

sturniraScatterplot <- 
ggplot(mergedWSymp, aes(y=FL, x=spbiog)) +
  geom_sina(aes(fill=species), pch=21, color="black", maxwidth=0.35) +
  xlab("")+ylab("forearm length")+
  stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, 
               geom="crossbar", width=0.4, size=0.4)+
  scale_x_discrete(labels=c(paste("S. hondurensis (A)", "\n", "n=", count(mergedWSymp,spbiog)[[2]][1]),
                            paste("S. hondurensis (S)", "\n", "n=", count(mergedWSymp,spbiog)[[2]][2]),
                            paste("S. parvidens (A)", "\n", "n=", count(mergedWSymp,spbiog)[[2]][3]),
                            paste("S. parvidens (S)", "\n", "n=", count(mergedWSymp,spbiog)[[2]][4])))+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
  geom_segment(aes(x = 1.5, y = 40, xend = 1.5, yend = 47),lineend = "round")+
  geom_segment(aes(x = 3.5, y = 37, xend = 3.5, yend = 45),lineend = "round")+
  guides(fill=FALSE)



########## FOR BEST

# Split by species and co occurrence 
hondSym <- mergedWSymp %>% filter(species=="hondurensis", sympatry=="both")
hondAll <- mergedWSymp %>% filter(species=="hondurensis", sympatry=="allop")
parvSym <- mergedWSymp %>% filter(species=="parvidens", sympatry=="both")
parvAll <- mergedWSymp %>% filter(species=="parvidens", sympatry=="allop")

# Run models for S. hondurensis (intraspecific)
BestHondSymAll <- BESTmcmc(hondSym$FL, hondAll$FL)
# optionally, save objects with MCMC results
saveRDS(object = BestHondSymAll, file="BESThondSymAll.rds")
# summarize results
summHondBEST <- summary(BestHondSymAll, credMass = 0.99) %>% as.data.frame()
# save plot object
hondDiffs <- plot(BestHondSymAll, col="#99C2CD", credMass = 0.99)
# wrangle object into something suitable for ggplot
hondDiffsDF <- data.frame(counts=hondDiffs$counts,
                          dens=hondDiffs$density,
                          mids=hondDiffs$mids)
# plot S. hondurensis
hondDiffsgrob <- 
ggplot(hondDiffsDF, aes(x=mids, y=dens))+geom_bar(stat="identity", fill="#335EAD")+
      geom_segment(aes(x = summHondBEST["muDiff", "HDIlo"],
                       y = 0, xend = summHondBEST["muDiff", "HDIup"], 
                       yend = 0), size=2.5)+
        xlab(expression(mu[1] - mu[2]))+ylab("density")+
        geom_segment(aes(x=0, y=0, xend=0, yend=max(hondDiffsDF$dens)*.75), linetype=2)+
        annotate("text", x=0, y=max(hondDiffsDF$dens)*.79,
                  label="0")
  

# Run models for S. parvidens (intraspecific)
BestParvSymAll <- BESTmcmc(parvSym$FL, parvAll$FL)
# optionally, save objects with MCMC results
saveRDS(object=BestParvSymAll, file="BestParvSymAll.rds")
# summarize results
summParvBEST <- summary(BestParvSymAll, credMass = 0.99) %>% as.data.frame()
# save plot object
parvDiffs <- plot(BestParvSymAll, col="#FFD271", credMass = 0.99)
# wrangle object into something suitable for ggplot
parvDiffsDF <- data.frame(counts=parvDiffs$counts,
                          dens=parvDiffs$density,
                          mids=parvDiffs$mids)

# plot S. parvidens
parvDiffsgrob <- 
  ggplot(parvDiffsDF, aes(x=mids, y=dens))+geom_bar(stat="identity", fill="#EEBB33")+
  geom_segment(aes(x = summParvBEST["muDiff", "HDIlo"],
                   y = 0, xend = summParvBEST["muDiff", "HDIup"], 
                   yend = 0), size=2.5)+
  xlab(expression(mu[1] - mu[2]))+ylab("density")+
  geom_segment(aes(x=0, y=0, xend=0, yend=max(parvDiffsDF$dens)*.75), linetype=2)+
  annotate("text", x=0, y=max(parvDiffsDF$dens)*.79,
           label="0")



# put plots together
grid.arrange(sturniraScatterplot, arrangeGrob(hondDiffsgrob, parvDiffsgrob, ncol=2), ncol=1, heights=c(1,0.52))
# save eps size was 820*700 px


####################################
# results table (intraspecific comparisons)
summHondBEST <- summary(BestHondSymAll, credMass = 0.99) %>% as.data.frame() %>% mutate(sp="hondurensis")
summParvBEST <- summary(BestParvSymAll, credMass = 0.99) %>% as.data.frame() %>% mutate(sp="parvidens")

#get the dimnames as row names for the df
dimnamesVec <- c(dimnames(summary(BestParvSymAll, credMass = 0.99) )[[1]])
sympAllBESToutput <- bind_rows(summHondBEST, summParvBEST)
sympAllBESToutput$param <- dimnamesVec #let it recycle
sympAllBESToutput <- sympAllBESToutput %>% select(param, everything())
#disk copy
write.csv(sympAllBESToutput, "LocdiffsBEST.csv")


#######  Overall differences 

# subset data
hondComp <- sturnFull %>% filter(species=="hondurensis")
parvComp <- sturnFull %>% filter(species=="parvidens")
# run model (interspecific)
BestOverallSp <- BESTmcmc(hondComp$FL, parvComp$FL)
# save object (optional)
saveRDS(object = BestOverallSp, file="bestOverall.rds")
plot(BestOverallSp, credMass = 0.99) #optional
# summarize
summALL <- summary(BestOverallSp, credMass = 0.99)%>% as.data.frame() %>% mutate(sp="both")
dimnamesVec <- c(dimnames(summary(BestOverallSp, credMass = 0.99) )[[1]])
summALL$param <- dimnamesVec 
summALL <- summALL %>% select(param, everything())

#disk copy
write.csv(summALL, "diffsALLsp.csv")

# further comparisons 
# run models Hondurensis vs Parvidens - sympatric (interspecific)
BestHondParvSym <- BESTmcmc(hondSym$FL, parvSym$FL)
round(summary(BestHondParvSym,credMass = 0.99),3)
# run models Hondurensis vs Parvidens - allopatric (interspecific)
BestHondParvAll <- BESTmcmc(hondAll$FL,parvAll$FL)
round(summary(BestHondParvAll, credMass = 0.99),3)


