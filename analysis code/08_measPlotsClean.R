# figure 2, descriptive plots
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggforce)
library(cowplot)
library(extrafont)
library(here)
library(rlang)

# read data
# fl and W
sturnfull <- read.csv(here("data","flwbiocl.csv"),stringsAsFactors = FALSE,encoding = "UTF-8")
# skulls
batskulls <- read.csv(here("data","cranBiocl.csv"),stringsAsFactors = FALSE,encoding = "UTF-8")

# cowplot settings
theme_set(theme_cowplot(font_size=12, font_family = "Roboto"))
loadfonts(device = "postscript")

# top panel, elevations
elevationsHist <- 
ggplot(sturnfull) +
  geom_density_ridges(stat="binline",aes(x=elevation, y=species,fill=species), color="black",
                      binwidth=100,scale=1) +  
    scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
    ggtitle("a)")+xlab("elevation (m)")+ ylab("frequency")+
    theme(plot.title =  element_text(hjust=0,size=12),
        text = element_text(size=12))+
  scale_y_discrete(expand = c(0, 0),labels=c("S. hondurensis","S. parvidens"))+
  scale_x_continuous(expand = c(0.01, 0)) +
  theme(legend.position="none",axis.text.y = element_text(face="italic"))

elevationsHist

# forearm length, overall
sturniraScatterplotFL <- 
  ggplot(sturnfull) +
  geom_sina(aes(y=FL, x=species,fill=species), pch=21, color="black", maxwidth=0.35, size=0.8) +  
  stat_summary(aes(y=FL, x=species),fun.y=mean, fun.ymin=mean, fun.ymax=mean, 
               geom="crossbar", width=0.4, size=0.3,color="black")+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
  ggtitle("b)")+xlab("")+ ylab("forearm length (mm)")+
  theme(plot.title =  element_text(hjust=0,size=12),
        text = element_text(size=12),
        axis.text.x = element_blank())+
  theme(legend.position="none")


######################################
# body mass, overall
sturniraScatterplotW <- 
  ggplot(sturnfull) +
  geom_sina(aes(y=W, x=species,fill=species), pch=21, color="black", maxwidth=0.35,size=0.8) +  
  stat_summary(aes(y=W, x=species),fun.y=mean, fun.ymin=mean, fun.ymax=mean, 
               geom="crossbar", width=0.4, size=0.3,color="black")+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
  ggtitle("c)")+xlab("")+ylab("body mass (g)")+
  theme(plot.title =  element_text(hjust=0,size=12),
        text = element_text(size=12),axis.text.x = element_blank())+
  theme(legend.position="none")


################################
#### skulls
# isosize fn
source(here("analysis code","Baur_isosize_fn.R"))

# skulls, but working only with complete cases
batskulls <- batskulls %>% dplyr::filter(complete.cases(.)) 

# measurements only
Skullsmat <-  batskulls %>% select(BRH:MSF)
#### calculate isosizes
batskulls$isosizes <- isosize(Skullsmat)[,1]

# isosize, overall
sturniraScatterplotIS <- 
  ggplot(batskulls) +
  geom_sina(aes(y=isosizes, x=species,fill=species), pch=21, color="black", maxwidth=0.35, size=0.8) +  
  stat_summary(aes(y=isosizes, x=species),fun.y=mean, fun.ymin=mean, fun.ymax=mean, 
               geom="crossbar", width=0.4, size=0.3,color="black")+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
  ggtitle("d)")+xlab("")+ylab("skull isosize")+
  theme(plot.title =  element_text(hjust=0,size=12),text = element_text(size=12),
        axis.text.x = element_blank())+
  theme(legend.position="none")#+

#########################
### head length
# read the data
batskulls <- read.csv(here("data","cranBiocl.csv"),stringsAsFactors = FALSE,encoding = "UTF-8") 

# head length, overall
sturniraScatterplotCBL <- 
  ggplot(batskulls) +
  geom_sina(aes(y= CBL, x=species,fill=species), pch=21, color="black", maxwidth=0.35,size=0.8) +  
  stat_summary(aes(y= CBL, x=species),fun.y=mean, fun.ymin=mean, fun.ymax=mean, 
               geom="crossbar", width=0.4, size=0.3,color="black")+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
  ggtitle("e)")+xlab("")+ylab("head length (mm)")+
  theme(plot.title =  element_text(hjust=0,size=12),text = element_text(size=12),
        axis.text.x = element_blank())+
  theme(legend.position="none")


######## side by side plots 
# dummy plot for shared legend
sturniraScatterplotLegend <- 
  ggplot(batskulls, aes(y=CBL, x=species)) +
  geom_sina(aes(fill=species),pch=21,color="black",maxwidth=0.35) +
  xlab("")+ylab("skull isosize")+
  scale_fill_manual(values = c("#335EAD","#EEBB33"),labels=c("S. hondurensis","S. parvidens"))+
  theme(legend.text = element_text(face = "italic"),legend.title=element_blank())

## extract the legend 
sharedLegend <- get_legend(sturniraScatterplotLegend+theme(legend.position="bottom",legend.justification="center"))


# fig 1
fig1 <- plot_grid(sturniraScatterplotFL,sturniraScatterplotW,sturniraScatterplotIS,sturniraScatterplotCBL, nrow = 1)
fig1

## Figure 1
# add the legend to the row we made earlier 
fig1L <- plot_grid(fig1, sharedLegend,ncol = 1, rel_widths = c(1,2),rel_heights = c(1.1, .08))
# draw
fig1L
# shared species 'label'
fig1Ls <- add_sub(fig1L, "species",vjust = -3.1,size=12)
# draw
fig1Ls <- ggdraw(fig1Ls)
fig1Ls

# fig, histograms at the top
fig1fin <- plot_grid(elevationsHist,fig1Ls,ncol=1,rel_heights = c(0.7, 1.3))
fig1fin

### export
ggsave(filename = here("figures","fig2.png"),fig1fin,width = 7.5, height=4, units = "in",dpi = 300)
