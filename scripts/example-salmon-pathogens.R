# Example of how to conduct a spatial analysis, using a dataset of pathogen 
# infections of Chinook salmon in British Columbia

# References: 

## 1. Load packages and data ----
library(tidyverse); theme_set(theme_bw())
library(sf)
library(readxl)
library(leaflet)
library(ggforce)
library(reshape2)

bc.coast <- read_sf("data/bc-coast.shp")

pathogens <- read_excel("data/disease-data.xlsx",
                              guess_max = 5719, sheet = "chinook")

## 2. Data exploration ----
glimpse(pathogens) # what variables does the dataset include, 5719 fish sampled

# 2A. Look at the Chinook salmon data

# plot the data spatially
ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = pathogens,
             aes(x = Longitude, y = Latitude), alpha = 0.3) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58))

# what is the temporal scope
table(pathogens$Year)

ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = subset(pathogens, Stock_Region == "ECVI"),
             aes(x = Longitude, y = Latitude), alpha = 0.3) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58)) +
  theme_bw() +
  facet_wrap(~ Year) # means "let's make a map for each year"

# different river populations were sampled. Let's take a look:
table(pathogens$Stock_Region) # check which rivers the Chinook came from
# (this was determined using a population genetics SNP chip)

ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = subset(pathogens, Stock_Region != "NA"), # exclude unknowns
             aes(x = Longitude, y = Latitude), alpha = 0.8) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58)) +
  facet_wrap(~ Stock_Region) # means let's make a map for each population

# Look at sample size distribution
# we will transform the coordinates from longitude/latitude to UTM zone 10
pathogens.wide.LL <- st_as_sf(pathogens, coords = c("Longitude", "Latitude"),
                              crs = st_crs(4326), remove = FALSE)
pathogens.wide.UTM <- st_transform(pathogens.wide.LL,
                                   crs = st_crs(32610)) # transform to UTM zone 10
st_crs(pathogens.wide.UTM)

# square grid with 30 km wide grid cells
grid <- st_make_grid(pathogens.wide.UTM, cellsize = 30000,
                     crs = st_crs(32610), what = "polygons")
bc.coast.UTM <- st_transform(bc.coast, crs = st_crs(32610))

ggplot() +
  geom_sf(data = bc.coast.UTM) +
  geom_sf(data = pathogens.wide.UTM) +
  geom_sf(data = grid, fill = NA) +
  coord_sf(datum = 4326, crs = 32610,
           xlim = c(-400000, 600000), ylim = c(5200000, 6600000))

# make a focused map of Vancouver Island
ggplot() +
  geom_sf(data = pathogens.wide.UTM, aes(colour = Ec_Wc_VanIsle)) 

van.island.UTM <- subset(pathogens.wide.UTM, Ec_Wc_VanIsle %in% c("ECVI", "WCVI") &
                                          Latitude < 51)
ggplot() +
  geom_sf(data = van.island.UTM, aes(colour = Ec_Wc_VanIsle)) 

# make a hexagonal grid # 10km grid spacing
hexagon.grid <- st_make_grid(van.island.UTM, cellsize = 10000,
                             crs = st_crs(32610), square = FALSE, what = "polygons")
class(hexagon.grid)
hexagon.grid <- st_as_sf(hexagon.grid)
hexagon.grid$grid_id <- seq(1, nrow(hexagon.grid))

ggplot() +
  geom_sf(data = van.island.UTM, aes(colour = Ec_Wc_VanIsle)) +
  geom_sf(data = hexagon.grid, fill = NA)

van.island.UTM <- st_join(x = van.island.UTM, y = hexagon.grid)

grid.sample.size <- data.frame(table(van.island.UTM$grid_id))
names(grid.sample.size) <- c("grid_id", "n")
grid.sample.size <- merge(x = hexagon.grid, grid.sample.size, by = "grid_id",
                           all.x = FALSE, all.y = FALSE)

ggplot() +
  geom_sf(data = bc.coast.UTM) +
  geom_sf(data = grid.sample.size, aes(fill = n)) +
  scale_fill_viridis_b(breaks = c(1, 10, 25, 50, 100, 802)) +
  coord_sf(datum = 4326, crs = 32610,
           xlim = c(50000, 530000), ylim = c(5300000, 5690000)) +
  labs(fill = "Sample\nsize")

## 2B. Look at the pathogens data ----

## what infectious agents were sampled
names(pathogens)[40:80] # 41, see table 1 in Bass et al. 2023 for more info

# spatial distribution of pathogens
# for ggplot, we need to convert the pathogens dataframe from "wide" to "long"

van.island.long <- pathogens %>% 
  subset(Ec_Wc_VanIsle %in% c("ECVI", "WCVI") & Latitude < 51) %>% 
  select(Fish, Species, Year, Stock_Region, Longitude, Latitude,
         arena1:vi_sal) %>%  # select key variables
  melt(., id.vars = c("Fish", "Species", "Year", "Stock_Region", "Longitude",
                      "Latitude"))
nrow(van.island.long) # 41 * 5060 = 207460, checks out!
names(van.island.long)[7:8] <- c("pathogen", "average_copy_number")

# since there are 41 pathogens, we will split this figure into two by using
# ggforce::facet_wrap_paginate()
str(van.island.long)

van.island.long$average_copy_number[van.island.long$average_copy_number %in% c("NA", "FALSE")] <- NA
van.island.long$average_copy_number <- as.numeric(van.island.long$average_copy_number)
van.island.long <- subset(van.island.long, !is.na(average_copy_number))

van.island.long$average_copy_number[van.island.long$average_copy_number < 0] <- 0
van.island.long.non.zero <- subset(van.island.long, average_copy_number > 0)

# plot distribution of pathogens
# using latitude/longitude (instead of UTM zone 10 coordiantes) since this is just
# a simple visualization
ggplot() +
  geom_sf(data = bc.coast) +
  # plot only non-zeros currently
  # log-transform to improve visualization
  geom_point(data = van.island.long.non.zero,
             aes(x = Longitude, y = Latitude, size = average_copy_number),
             alpha = 0.8) +
  scale_size_continuous(breaks = c(1, 1000, 1000000, 1e+8)) +
  coord_sf(xlim = c(-130.5, -122), ylim = c(47.6, 51.4)) +
  facet_wrap_paginate(~ pathogen, nrow = 3, ncol = 3, page = 1) +
  theme(legend.position = "bottom")
# change page to look at all the different pathogens

# 2C: evaluate regional pathogen load ----
# load Fisheries and Oceans Canada Pacific Fishery Management Areas (PFMAs)
# I saved shapefile as R Data Storage file to reduce file size for workshop
sf_use_s2(FALSE)
pfma <- readRDS("data/DFO-management-areas.RDS") 

table(sf::st_is_valid(pfma)) # this

pathogens.sf <- st_as_sf(van.island.long, coords = c("Longitude", "Latitude"),
                         crs = st_crs(4326), remove = FALSE)
# CRS 4326 is the coordinate reference system for WGS 84 Mercator

st_crs(pfma) # crs 4326
st_crs(pathogens.sf) # crs 3857

pfma.sf <- st_transform(pfma, crs = st_crs(4326))
sf_use_s2(FALSE)
pathogens.pfma <- pathogens.sf %>% 
  st_join(x = ., y = pfma.sf)

# note that some samples did not fall within in a PFMA
nrow(subset(pathogens.pfma, is.na(MGNT_AREA))) # 9702
nrow(subset(pathogens, is.na(Longitude)))

# these samples are mostly in US waters. Samples from Canada that failed to
# join with the PFMA layer are likely GPS errors
no.pfma <- subset(pathogens.pfma, is.na(MGNT_AREA))

ggplot() +
  geom_sf(data = pfma.sf, fill = "lightblue", colour = "lightblue") +
  geom_sf(data = bc.coast) + 
  geom_point(data = no.pfma,
             aes(x = Longitude, y = Latitude, colour = MGNT_AREA)) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58))

ggplot() +
  geom_sf(data = pfma.sf, fill = "lightblue", colour = "lightblue") +
  geom_sf(data = bc.coast) + 
  geom_point(data = no.pfma,
             aes(x = Longitude, y = Latitude, colour = MGNT_AREA)) +
  coord_sf(xlim = c(-125.5, -125), ylim = c(50.2, 50.5)) 

# if you were working up an analysis, you would make sure to correct the coordinates
# Since this is just an example, we will remove samples that are on land, as well
# as those from American waters

# assess regional pathogen load
names(pathogens.pfma)
pathogens.pfma <- pathogens.pfma %>% 
  # remove samples from US or faulty coordinates (i.e. on land):
  subset(!is.na(MGNT_AREA)) %>% 
  # make region variable: note \n indicates a line break
  mutate(region = case_when(MGNT_AREA %in% c(1:5, 101:104, 108, 109, 142) ~ "Haida/\nNorth Coast", # Haida Gwaii/North  Coast
                            MGNT_AREA %in% c(6:10, 11, 106, 107, 110) ~ "Central\nCoast",
                            MGNT_AREA %in% c(12:19, 28, 29) ~ "East Van. I.",
                            MGNT_AREA %in% c(20:27, 111, 121:127) ~ "West Van. I.")) %>% 
  mutate
table(pathogens.pfma$pathogen)
ic <- subset(pathogens.pfma, pathogen == "ic_hof")

ggplot() +
  geom_sf(data = bc.coast) + 
  geom_point(data = ic,
             aes(x = Longitude, y = Latitude, colour = region)) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58))

ggplot() +
  geom_boxplot(data = pathogens.pfma,
               aes(x = region, y = average_copy_number, fill = region),
               show.legend = FALSE) + 
  facet_wrap_paginate(~ pathogen, nrow = 3, ncol = 3, page = 1, scales = "free") +
  labs(x = "Region", y = "Average copy number")
# change page to move through different options

# using a log scale can help interpretation

ggplot() +
  geom_boxplot(data = pathogens.pfma,
               aes(x = region, y = log(average_copy_number + 1), fill = region),
               show.legend = FALSE) + 
  facet_wrap_paginate(~ pathogen, nrow = 3, ncol = 3, page = 3, scales = "free") +
  labs(x = "Region", y = "Log average copy number + 1")

## 3. Spatial model ----


## make an interactive map
leaflet(pathogens) %>% addTiles() %>% 
  addMarkers(lng = ~Longitude, lat = ~Latitude,
             label = ~Stock)
