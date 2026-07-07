##########
#Figure 1
##########
library(dplyr)
library(ggmap)
library(scatterpie)
library(wesanderson)
library(PNWColors)
library(MetBrewer)
library(maps)
library(mapdata)
#library(maptools)  #for shapefiles
library(scales)  #for transparency
#library(geodata)

#cg<-read.table("~/3waymerged_sampleinfo.txt",sep="\t",na.strings = c("","NA"),header=T)
cg <- read.table("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/Chp1_DroughtandCommonGardenInfo_Merged-MetaData.tsv", sep="\t",na.strings = c("","NA"),header=T)
head(cg)

#assign a dataset variable
cg$dataset<-"Chp1"
cg$dataset[is.na(cg$Sample_Name_Notes) == T]<-"CommonGarden"
cg$dataset[is.na(cg$Sample_Name_Notes) == F]<-"Drought"


#get pop level coords
#commented bc I didnt have state var and I dont want to add it
#coord_byenv<-cg %>% group_by(Environment, Long, Lat, dataset) %>% dplyr:::summarise(n=n(), state=first(state))
coord_byenv<-cg %>% group_by(Environment, Long, Lat, dataset) %>% dplyr:::summarise(n=n())

#get map base
library(raster)
#states    <- c("New York","Pennsylvania","Maryland","West Virginia","Virginia","Kentucky","Ohio", "Michigan","Indiana","Illinois","Wisconsin","Minnesota","Iowa","Missouri","Kansas","Nebraska","South Dakota","North Carolina","Tennessee","Mississippi", "Oklahoma","Lake Michigan","Lake Ontario","Lake Superior")

states    <- c("Kentucky","Ohio","Indiana","Illinois","Iowa","Missouri","Kansas","Nebraska","South Dakota","North Carolina","Tennessee","Mississippi", "Oklahoma","Lake Michigan","Lake Ontario","Lake Superior")

provinces <- c("Ontario")

us <- gadm(country="USA",level=1)
canada <- gadm(country="CAN",level=1)

#us.states <- us[us$NAME_1 %in% states,]
#ca.provinces <- canada[canada$NAME_1 %in% provinces,]

#last two where giving plotting error
us.states <- map_data("state")
ca.provinces <- map_data("world", region = "Canada")

#lakes
lakes <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'lakes', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

#ocean
ocean <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'ocean', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

# rivers
rivers <- rnaturalearth::ne_download(scale = 10, 
                                     type = 'rivers_lake_centerlines', 
                                     category = 'physical')  %>% 
  sf::st_as_sf(lakes110, crs = 4269)

library("rnaturalearth")
#install.packages("rnaturalearthhires")
library("rnaturalearthhires")
library(rnaturalearth)
library(dplyr)
library(raster)
library(sf)
library(tidyverse)
library(ggrepel)

#get US states outlines
in_sf <- ne_states(geounit = "United States of America",
                   returnclass = "sf")

uss <- in_sf %>%
  mutate(
    lon = st_coordinates(st_centroid(geometry))[, 1],
    lat  = st_coordinates(st_centroid(geometry))[, 2]
  )

#move Michigan label and add Ontario label
#uss[uss$NAME_1=="Michigan",]$lon<--85.0554
#uss[uss$NAME_1=="Michigan",]$lat<-44.00902
#uss<-uss %>% add_row(NAME_1 = "Ontario", lon = -78.5554, lat=45)

#rename environmental variables for legend
coord_byenv$Environment[coord_byenv$Environment == 2] <- "Agricultural"
coord_byenv$Environment[coord_byenv$Environment == 1] <- "Natural"
cg$env[cg$Environment == 2] <- "Agricultural"
cg$env[cg$Environment == 1] <- "Natural"

#cg<-cg[complete.cases(cg$env),]

coord_byenv$Environment<-as.factor(coord_byenv$Environment)
cg$Environment<-as.factor(cg$Environment)
cg$Environment<-factor(cg$Environment, levels = c("Agricultural", "Natural"))


#generate full map base
library(ggthemes)
plain1<- 
  ggplot()+
  geom_path(data=us.states,aes(x=long,y=lat,group=group))+
  geom_path(data=ca.provinces, aes(x=long,y=lat,group=group))+
  coord_map() +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          color = "black")  +
  geom_sf(data = ocean,
          mapping = aes(geometry = geometry),
          color = "black")  +
  geom_sf(data = rivers,
          mapping = aes(geometry = geometry),
          color = "grey80",alpha=.75)  +
  #theme_nothing() +
  scale_x_continuous(limits=c(-97,-83)) +
  scale_y_continuous(limits=c(38,43)) +
  coord_sf(xlim=c(-97,-83)) +
  geom_text(
    data = uss,
    aes(x = lon,
        y = lat,
        label = name),cex=3, alpha=.8) 


#import color palettes
bay<-pnw_palette("Bay",7,type="continuous")
gaug<-met.brewer(name="Gauguin")

#package for multiple color schemes
library(ggnewscale)

#rename environmental variables for legend
coord_byenv$Environment[coord_byenv$Environment == 2] <- "Agricultural"
coord_byenv$Environment[coord_byenv$Environment == 1] <- "Natural"
cg<-cg[complete.cases(cg$Environment),]

#plot map base with data points
#change color back by commenting out scale_color_manual
final <- plain1 + 
  geom_jitter(data=coord_byenv,aes(y = Lat, x = Long, color= Environment, size=n, shape=Environment), alpha=.6, height=.2, width=.15) +
  #scale_color_manual(values=c("#60CEACFF","#382A54FF")) +
  scale_color_manual(values=c("gray20","gray60")) +
  scale_shape_manual(values=c("circle","triangle")) +
  #new_scale_color() +
  #geom_jitter(data=cg[cg$dataset == "Herbarium",],aes(y = lat, x = long, color=year, shape=env), size=3, height=.1,width=.05, alpha=.9) +
  #scale_color_viridis_c() +
  xlab("Longitude") +
  ylab("Latitude")

ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Presentations/Figures/Chp1_map.png", plot = final, device = "png")




##########
#Range Maps
##########

library(data.table)
library(tidyverse)
#lat: (29,42)
#long: c(-95,-125)) +
palmeri <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/palmeri.csv")
palmeri <- palmeri %>%
  filter(countryCode == "US") %>%
  #filter(decimalLatitude > 29 & decimalLatitude < 42 &  decimalLongitude > -125 & decimalLongitude < -95) %>% 
  select(species, decimalLatitude, decimalLongitude)

#palmeri_outlier <- filter(palmeri, decimalLongitude > -115 & decimalLongitude < -108 & decimalLatitude > 37)
#remove utah outlier 40.76078, -111.8910
#palmeri <- palmeri[-2254,]

albus <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/albus.csv")
albus <- albus %>%
  filter(countryCode == "US") %>%
  #filter(decimalLatitude > 29 & decimalLatitude < 42 &  decimalLongitude > -125 & decimalLongitude < -95) %>%
  select(species, decimalLatitude, decimalLongitude)

retroflexus <- fread("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Coding/Rotation_Winter2025/Data/retroflexus.csv")
retroflexus <- retroflexus %>%
  filter(countryCode == "US") %>%
  #filter(decimalLatitude > 29 & decimalLatitude < 42 &  decimalLongitude > -125 & decimalLongitude < -95) %>%
  select(species, decimalLatitude, decimalLongitude)

#get map base
library(raster)
library(geodata)

states<- c("California", "Arizona", "New Mexico", "Texas")

us <- gadm(country="USA",level=1)

#last two where giving plotting error
us.states <- map_data("state")

#lakes
lakes <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'lakes', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

#ocean
ocean <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'ocean', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

# rivers
rivers <- rnaturalearth::ne_download(scale = 10, 
                                     type = 'rivers_lake_centerlines', 
                                     category = 'physical')  %>% 
  sf::st_as_sf(lakes110, crs = 4269)

library("rnaturalearth")
library("rnaturalearthhires")
library(rnaturalearth)
library(sf)
library(ggrepel)

#get US states outlines
in_sf <- ne_states(geounit = "United States of America",
                   returnclass = "sf")

uss <- in_sf %>%
  mutate(
    lon = st_coordinates(st_centroid(geometry))[, 1],
    lat  = st_coordinates(st_centroid(geometry))[, 2]
  )


#generate full map base
library(ggthemes)
plain1<- 
  ggplot()+
  geom_path(data=us.states,aes(x=long,y=lat,group=group))+
  coord_map() +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          color = "black")  +
  geom_sf(data = ocean,
          mapping = aes(geometry = geometry),
          color = "black")  +
  geom_sf(data = rivers,
          mapping = aes(geometry = geometry),
          color = "grey80",alpha=.75)  +
  #theme_nothing() +
  scale_x_continuous(limits=c(-95,-125)) +
  scale_y_continuous(limits=c(29,42)) +
  coord_sf(xlim=c(-95,-125)) +
  geom_text(
    data = uss,
    aes(x = lon,
        y = lat,
        label = name),cex=3, alpha=.8) 


#import color palettes
library(PNWColors)
library(MetBrewer)
bay<-pnw_palette("Bay",7,type="continuous")
gaug<-met.brewer(name="Gauguin")

#create range polygon
palmeri_points_sf <- st_as_sf(palmeri, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
# Convex Hull method (smallest convex polygon containing points)
boundary_palmeri <- st_convex_hull(st_union(palmeri_points_sf ))
# Get points that intersect the boundary polygon
#boundary_points <- points_sf[st_intersects(boundary, points_sf)[[1]], ]
albus_points_sf <- st_as_sf(albus, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
boundary_albus <- st_convex_hull(st_union(albus_points_sf))
retro_points_sf <- st_as_sf(retroflexus, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
boundary_retro <- st_convex_hull(st_union(retro_points_sf))


#plot map base with data points
#change color back by commenting out scale_color_manual
#colors <- c("#0B0405FF", "#382A54FF", "#395D9CFF", "#3497A9FF", "#60CEACFF", "#DEF5E5FF")
#final_range <- plain1 + geom_sf(data = boundary_palmeri, aes(fill = "A. palmeri"), color = "#382A54FF", size = 1, alpha = .2) + geom_sf(data = boundary_albus, aes(fill = "A. albus"), color = "#395D9CFF", size = 1, alpha = .2) + geom_sf(data = boundary_retro, aes(fill = "A. retroflexus"), color = "#3497A9FF", size = 1, alpha = .2) +
  xlab("Longitude") +
  ylab("Latitude") + scale_fill_manual(name = "Species", values=c("A. palmeri" = "#382A54FF", "A. albus"= "#395D9CFF", "A. retroflexus" = "#3497A9FF")) 

# Background points geom_sf(data = boundary_palmeri, fill = NA, color = "red", size = 1)

  final_range <-  plain1 + geom_jitter(data=palmeri, aes(y = decimalLatitude, x = decimalLongitude, color= species), size = .8, alpha=.3, height=.1, width=.15) + geom_jitter(data=albus,aes(y = decimalLatitude, x = decimalLongitude, color= species), size = .8, alpha=.3, height=.1, width=.15) + geom_jitter(data=retroflexus,aes(y = decimalLatitude, x = decimalLongitude, color= species), size = .8, alpha=.3, height=.1, width=.15) + xlab("Longitude") + ylab("Latitude") + scale_color_manual(name = "Species", values=c("#382A54FF", "#395D9CFF", "#3497A9FF"))

ggsave("/Users/libbypolston/Desktop/UChicago/Kreiner_lab/Presentations/Figures/Chp2_rangemaps.png", plot = final_range, device = "png")

#used this to get darker legend colors
plain1 + geom_jitter(data=palmeri, aes(y = decimalLatitude, x = decimalLongitude, color= species), size = .8, alpha=.8, height=.1, width=.15) + geom_jitter(data=albus,aes(y = decimalLatitude, x = decimalLongitude, color= species), size = .8, alpha=.8, height=.1, width=.15) + geom_jitter(data=retroflexus,aes(y = decimalLatitude, x = decimalLongitude, color= species), size = .8, alpha=.8, height=.1, width=.15) + xlab("Longitude") + ylab("Latitude") + scale_color_manual(name = "Species", values=c("#382A54FF", "#395D9CFF", "#3497A9FF"))
