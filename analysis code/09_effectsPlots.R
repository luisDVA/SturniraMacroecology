# sympatry effects plots
library(dplyr)
library(purrr)
library(here)
library(ggplot2)
library(cowplot)
library(extrafont)


# read the exported effects for plotting
Alleffects <- 
list.files(here("lmmEffects"), pattern = ".csv$", full.names = TRUE) %>% 
  map_df(read.csv)

# cowplot settings
theme_set(theme_cowplot(font_family = "Roboto",font_size = 12))

# object for each measurement  
# forearm
fleff <- 
Alleffects %>% filter(measurement=="forearm") %>%
ggplot(aes(x=species, y=fit, pch=species,fill=species),color="black") + 
  geom_pointrange(aes(ymin=lower, ymax=upper),size=0.9)+
  facet_wrap(~sympatry,labeller=labeller(sympatry = c(allop="allopatry",both="sympatry")))+
  ylab("forearm length (mm)")+xlab("")+
  scale_shape_manual(values = c(21, 23))+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
  theme(legend.position="none",axis.text.x = element_blank(),
        strip.background = element_rect(fill ="#FDFCFF"))

# body mass
bmeff <- 
Alleffects %>% filter(measurement=="bodyMass") %>%
  ggplot(aes(x=species, y=fit, pch=species,fill=species),color="black") + 
  geom_pointrange(aes(ymin=lower, ymax=upper),size=0.9)+
  facet_wrap(~sympatry,labeller=labeller(sympatry = c(allop="allopatry",both="sympatry")))+
  ylab("body mass (g)")+xlab("")+
  scale_shape_manual(values = c(21, 23))+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
  theme(legend.position="none",axis.text.x = element_blank(),
        strip.background = element_rect(fill ="#FDFCFF"))

# skull isosize
ISeff <- 
Alleffects %>% filter(measurement=="isosize") %>%
  ggplot(aes(x=species, y=fit, pch=species,fill=species),color="black") + 
  geom_pointrange(aes(ymin=lower, ymax=upper),size=0.9)+
  facet_wrap(~sympatry,labeller=labeller(sympatry = c(allop="allopatry",both="sympatry")))+
  ylab("skull isosize")+xlab("")+
  scale_shape_manual(values = c(21, 23))+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
  theme(legend.position="none",axis.text.x = element_blank(),
        strip.background = element_rect(fill ="#FDFCFF"))

# head length
CBLeff <- 
Alleffects %>% filter(measurement=="CBL") %>%
  ggplot(aes(x=species, y=fit, pch=species,fill=species),color="black") + 
  geom_pointrange(aes(ymin=lower, ymax=upper),size=0.9)+
  facet_wrap(~sympatry,labeller=labeller(sympatry = c(allop="allopatry",both="sympatry")))+
  ylab("head length (mm)")+xlab("")+
  scale_shape_manual(values = c(21, 23))+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"))+
  theme(legend.position="none",axis.text.x = element_blank(),
        strip.background = element_rect(fill ="#FDFCFF"))


######## side by side plots 
# dummy plot for shared legend
# hack to separate the elements
legplot <- 
  Alleffects %>% filter(measurement=="CBL") %>%
  ggplot(aes(x=species, y=fit, pch=species,fill=species),color="black") + 
  geom_pointrange(aes(ymin=lower, ymax=upper),size=0.9)+
  facet_wrap(~sympatry,labeller=labeller(sympatry = c(allop="allopatry",both="sympatry")))+
  scale_shape_manual(values = c(21, 23),labels=c("S. hondurensis   ","S. parvidens"))+
  scale_fill_manual(values = c("#335EAD", "#EEBB33"),labels=c("S. hondurensis   ","S. parvidens"))+
  theme(legend.text = element_text(face = "italic"),legend.title=element_blank())

## extract the legend 
sharedLegend <- get_legend(legplot+theme(legend.position="bottom",legend.justification="center"))

# fig 3
fig3 <- plot_grid(fleff,bmeff,ISeff,CBLeff,nrow = 1)
fig3

## Figure 3
# add the legend to the row we made earlier. Give it one-third of the width of one plot (via rel_widths).
fig3L <- plot_grid(fig3, sharedLegend,ncol = 1, rel_heights = c(1.1, .08))
# draw
fig3L
# shared species 'label'
fig3Ls <- add_sub(fig3L, "species",vjust = -4.2,size=12)
# draw
fig3Ls <- ggdraw(fig3Ls)
fig3Ls

#export
ggsave(fig3Ls,filename = here("figures","fig3.png"), width = 8, height = 3.6, units = "in", dpi = 300)

