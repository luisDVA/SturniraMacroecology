# figure 1, map
library(sf)
library(rgdal)
library(dplyr)
library(here)
library(rnaturalearth)
library(tmaptools)
library(SDMTools)
library(ggplot2) # dev version that has geom_sf 

## read coordinates
# forearms
sturnFull <- read.csv(here("data","flwbiocl.csv"),stringsAsFactors = FALSE)
sturnpoints <- sturnFull %>% dplyr::select(species,longitude,latitude,FL)
coordinates(sturnpoints) <- ~longitude + latitude
proj4string(sturnpoints) <- CRS("+proj=longlat +datum=WGS84")
# skulls
batSkulls <- read.csv(here("data","cranBiocl.csv"), stringsAsFactors = FALSE)
skullspoints <- batSkulls %>% dplyr::select(species,longitude,latitude,CBL)
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

# copy of points in projected system (INEGI Lambert conformal)
sturnpointsproj <- spTransform(sturnpoints,CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +a=6378137 +b=6378136.027241431 +units=m +no_defs"))
ptdens <- smooth_map(sturnpointsproj)
ptdensities <- extract.data(sturnpointsproj@coords,ptdens$raster)

skullspointsproj <- spTransform(skullspoints,CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +a=6378137 +b=6378136.027241431 +units=m +no_defs"))
skullsptdens <- smooth_map(skullspointsproj)
skptdensities <- extract.data(skullspointsproj@coords,skullsptdens$raster)

# bind with points sf dataframe
sturnpointssf <-  cbind(sturnpointssf,ptdensities)
skullspointssf <- cbind(skullspointssf, skptdensities)

## plotting

# rescale sizes
sturnpointssf$FLsc <- round(scales::rescale(sturnpointssf$FL,c(1.5,8)),1)
skullspointssf$CBLsc <- round(scales::rescale(skullspointssf$CBL,c(1.5,8)),1)

## plot

fig1 <- 
ggplot()+
  geom_sf(data=divpolsf,fill= NA)+
  geom_sf(data=sturnpointssf,aes(shape=species,color="charcoal", fill=species,size=FLsc, alpha=1/ptdensities),
          show.legend = FALSE)+
  scale_shape_manual(values=c(21,23))+
  scale_fill_manual(values=c("#335EAD","#EEBB33"))+
  xlim(c(-108,-87))+ylim(c(14,24))+
  scale_alpha(range = c(.09, .8))+
  scale_size_identity()+
  geom_sf(data=skullspointssf,pch=4,color="black",aes(size=CBLsc,alpha=1/skptdensities),
          show.legend = FALSE)+
  theme(panel.grid.major = element_line(colour = 'transparent'),
        rect = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))




# export
ggsave(fig1,filename = here("figures","fig1.png"),width = 7, height = 4,units = "in",dpi=300)

