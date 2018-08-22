# figure 1, map
library(sf)
library(rgdal)
library(dplyr)
library(here)
library(rnaturalearth)
library(pointdensityP)
library(ggplot2) 

## read coordinates
# forearms
sturnFull <- read.csv(here("data","flwbiocl.csv"),stringsAsFactors = FALSE)
sturnpoints <- sturnFull %>% dplyr::select(species,longitude,latitude,FL)
# densities
sturnptdens <- pointdensity(sturnpoints,lat_col = "latitude",
                         lon_col = "longitude", grid_size = 40,radius = 80)
# merge
sturnpoints <- bind_cols(sturnpoints,sturnptdens) %>% rename(ptdensities=count)
sturnpoints <- sturnpoints %>% mutate(ptdensitiesSc=scales::rescale(ptdensities,c(0.0001,1)))
# set up sp
coordinates(sturnpoints) <- ~longitude + latitude
proj4string(sturnpoints) <- CRS("+proj=longlat +datum=WGS84")
# skulls
batSkulls <- read.csv(here("data","cranBiocl.csv"), stringsAsFactors = FALSE)
skullspoints <- batSkulls %>% dplyr::select(species,longitude,latitude,CBL)
skullptdens <- pointdensity(skullspoints,lat_col = "latitude",
                            lon_col = "longitude", grid_size = 40,radius = 80)
# merge
skullspoints <- bind_cols(skullspoints,skullptdens) %>% rename(ptdensities=count)
skullspoints <- skullspoints %>% mutate(ptdensitiesSc=scales::rescale(ptdensities,c(0.0001,1)))
# set up sp
coordinates(skullspoints) <- ~longitude + latitude
proj4string(skullspoints) <- CRS("+proj=longlat +datum=WGS84")
# convert to simple features
sturnpointssf <- st_as_sf(sturnpoints)
skullspointssf <- st_as_sf(skullspoints)

# basemap
divpol <- ne_countries(country = c("mexico","guatemala",
                                   "belize","el salvador",
                                   "honduras"),scale = 50)
divpolsf <- st_as_sf(divpol)


## plotting

# rescale sizes
sturnpointssf$FLsc <- round(scales::rescale(sturnpointssf$FL,c(1.5,8)),1)
skullspointssf$CBLsc <- round(scales::rescale(skullspointssf$CBL,c(1.5,8)),1)
## plot
fig1 <- 
  ggplot()+
  geom_sf(data=divpolsf,fill= NA)+
  geom_sf(data=sturnpointssf,aes(shape=species,color="charcoal", 
                                 fill=species,size=FLsc, alpha=ptdensitiesSc),
          show.legend = FALSE)+
  scale_shape_manual(values=c(21,23))+
  scale_fill_manual(values=c("#335EAD","#EEBB33"))+
  xlim(c(-108,-87))+ylim(c(14,24))+
  scale_alpha(range = c(.09, .8))+
  scale_size_identity()+
  geom_sf(data=skullspointssf,pch=4,color="black",aes(size=CBLsc,alpha=ptdensitiesSc),
          show.legend = FALSE)+
  theme(panel.grid.major = element_line(colour = 'transparent'),
        rect = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# export
ggsave(fig1,filename = here("figures","fig1hires.png"),width = 7, height = 4,units = "in",dpi=900)


