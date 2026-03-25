

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

# should make a grid for a sample size map

## 2B. Look at the pathogens data ----

## what infectious agents were sampled
names(pathogens)[40:80] # 41, see table 1 in Bass et al. 2023 for more info

# spatial distribution of pathogens
# for ggplot, we need to convert the pathogens dataframe from "wide" to "long"

pathogens.long <- pathogens %>% 
  select(Fish, Species, Year, Stock_Region, Longitude, Latitude,
         arena1:vi_sal) %>%  # select key variables
  melt(., id.vars = c("Fish", "Species", "Year", "Stock_Region", "Longitude",
                      "Latitude"))
nrow(pathogens.long) # 41 * 5719 = 234479, checks out!
names(pathogens.long)[7:8] <- c("pathogen", "average_copy_number")

# since there are 41 pathogens, we will split this figure into two by using
# ggforce::facet_wrap_paginate()
str(pathogens.long)

pathogens.long$average_copy_number[pathogens.long$average_copy_number %in% c("NA", "FALSE")] <- NA
pathogens.long$average_copy_number <- as.numeric(pathogens.long$average_copy_number)
pathogens.long <- subset(pathogens.long, !is.na(average_copy_number))

pathogens.long$average_copy_number[pathogens.long$average_copy_number < 0] <- 0

pathogens.long.non.zero <- subset(pathogens.long, average_copy_number > 0)

ggplot() +
  geom_sf(data = bc.coast) +
  # plot only non-zeros currently
  geom_point(data = pathogens.long.non.zero,
             aes(x = Longitude, y = Latitude, size = average_copy_number),
             alpha = 0.8) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58)) +
  facet_wrap_paginate(~ pathogen, nrow = 3, ncol = 3, page = 1) +
  theme(legend.position = "bottom")
# change page to look at all the different pathogens

# 2C: evaluate regional pathogen load ----
# load Fisheries and Oceans Canada Pacific Fishery Management Areas (PFMAs)
# I saved shapefile as R Data Storage file to reduce file size for workshop
sf_use_s2(FALSE)
pfma <- readRDS("data/DFO-management-areas.RDS") 

table(sf::st_is_valid(pfma)) # this

pathogens.sf <- st_as_sf(pathogens.long, coords = c("Longitude", "Latitude"),
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
