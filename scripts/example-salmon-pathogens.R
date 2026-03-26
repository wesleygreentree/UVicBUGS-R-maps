# Example of how to conduct a spatial analysis, using a dataset of pathogen 
# infections of Chinook salmon in British Columbia

# References: 
# Dataset from:
# Bass et al. 2023. The spatial distribution of infectious agents in wild Pacific
# salmon along the British Columbia coast. Scientific Reports. 
# https://doi.org/10.1038%2Fs41598-023-32583-8

# Replicating analyses from:
# Bass et al. 2024. Intrinsic and extrinsic factors associated with the spatio-
# temporal distribution of infectious agents in early marine Chinook and coho
# salmon. Marine Ecology Progress Series. https://doi.org/10.3354/meps14581


## make an interactive map
leaflet(pathogens) %>% addTiles() %>% 
  addMarkers(lng = ~Longitude, lat = ~Latitude,
             label = ~Stock)

## 1. Load packages and data ----
library(tidyverse); theme_set(theme_bw())
library(sf)
library(readxl)
library(leaflet)
library(ggforce)
library(reshape2)
library(sdmTMB)

bc.coast <- read_sf("data/bc-coast.shp")

pathogens <- read_excel("data/disease-data.xlsx",
                              guess_max = 5719, sheet = "chinook")

## 2. Data exploration ----
glimpse(pathogens) # what variables does the dataset include, 5719 fish sampled

# 2A. Look at the Chinook salmon data ----

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
  geom_point(data = pathogens,
             aes(x = Longitude, y = Latitude), alpha = 0.8) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58)) +
  theme_bw() +
  facet_wrap(~ Year) # means "let's make a map for each year"

# different river populations were sampled. Let's take a look:
table(pathogens$Stock_Region) # check which rivers the Chinook came from
# (this was determined using a population genetics SNP chip)
table(pathogens$CU_name) # DFO conservation units, based on genetics

ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = subset(pathogens, Stock_Region != "NA"), # exclude unknowns
             aes(x = Longitude, y = Latitude), alpha = 0.8) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58)) +
  facet_wrap(~ Stock_Region) # means let's make a map for each population

# look at one example of differences in stock distribution: Columbia R and
# East Coast Vancouver Island populations have nearly opposite distribution
ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = subset(pathogens, Stock_Region %in% c("Columbia", "ECVI")), 
             aes(x = Longitude, y = Latitude), alpha = 0.8) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58)) +
  facet_wrap(~ Stock_Region) # means let's make a map for each population

# Look at sample size distribution
# we will transform the coordinates from longitude/latitude to UTM zone 10

# first need to make an sf object for pathogens using longitude/latitude
pathogens.LL <- st_as_sf(pathogens, coords = c("Longitude", "Latitude"),
                              crs = st_crs(4326), remove = FALSE)
pathogens.UTM <- st_transform(pathogens.LL,
                                   crs = st_crs(32610)) # transform to UTM zone 10
st_crs(pathogens.UTM) # note that units are metres

# square grid with 30 km wide grid cells
grid <- st_make_grid(pathogens.UTM, cellsize = 30000,
                     crs = st_crs(32610), what = "polygons", square = TRUE)
bc.coast.UTM <- st_transform(bc.coast, crs = st_crs(32610))

# take a look at the grid
ggplot() +
  geom_sf(data = bc.coast.UTM) +
  geom_sf(data = pathogens.UTM) +
  geom_sf(data = grid, fill = NA) +
  coord_sf(datum = 4326, crs = 32610,
           xlim = c(-400000, 600000), ylim = c(5200000, 6600000))

# make a focused map of Vancouver Island
ggplot() +
  geom_sf(data = pathogens.UTM, aes(colour = Ec_Wc_VanIsle)) 
van.island.UTM <- subset(pathogens.UTM, Ec_Wc_VanIsle %in% c("ECVI", "WCVI") &
                                          Latitude < 51)
ggplot() +
  geom_sf(data = van.island.UTM, aes(colour = Ec_Wc_VanIsle)) 

# make a hexagonal grid # 10km grid spacing
# replicates Figure 1 from Bass et al. 2024
hexagon.grid <- st_make_grid(van.island.UTM, cellsize = 10000,
                             crs = st_crs(32610), square = FALSE, what = "polygons")
# square = FALSE: hexagons; square = TRUE: square grid cells

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
str(van.island.long) # need to make average_copy_number a numeric variable

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
  scale_x_continuous(breaks = c(-130, -126, -123)) +
  scale_y_continuous(breaks = c(48, 49, 51)) +
  facet_wrap_paginate(~ pathogen, nrow = 3, ncol = 3, page = 1) +
  theme(legend.position = "bottom")
# change page to look at all the different pathogens

## 2C. Calculate prevalence in each grid cell with at least ----

lo.sal <- data.frame(matrix(data = NA, nrow = 0, ncol = 3))
names(lo.sal) <- c("grid_id", "total_n", "prevalence")

grid.vector <- sort(unique(van.island.UTM$grid_id))

for (i in grid.vector) {
  
  grid.id <- i
  
  grid.subset <- subset(van.island.UTM, grid_id == i)
  
  n <- nrow(grid.subset)
  positive <- nrow(subset(grid.subset, lo_sal > 0))  
  prop <- positive / n

  out <- data.frame(grid_id = grid.id,
                    total_n = n, 
                    prevalence = prop)
  
  lo.sal <- rbind(lo.sal, out)  
}

length(grid.vector)
nrow(lo.sal) # 268

lo.sal <- merge(x = hexagon.grid, y = lo.sal, by = "grid_id",
                          all.x = FALSE, all.y = FALSE)
lo.sal$Percentage <- lo.sal$prevalence * 100

ggplot() +
  geom_sf(data = bc.coast.UTM) +
  geom_sf(data = lo.sal, aes(fill = Percentage)) +
  scale_fill_viridis_c(option = "B", limits = c(0, 100)) +
  coord_sf(datum = 4326, crs = 32610,
           xlim = c(50000, 530000), ylim = c(5300000, 5690000))

# plot only grid cells with at least 10 samples
ggplot() +
  geom_sf(data = bc.coast.UTM) +
  geom_sf(data = subset(lo.sal, total_n >= 5), aes(fill = Percentage)) +
  scale_fill_viridis_c(option = "B", limits = c(0, 100)) +
  coord_sf(datum = 4326, crs = 32610,
           xlim = c(50000, 530000), ylim = c(5300000, 5690000))
ggsave("figures/dev/loma-salmonae-n5.PNG",
       width = 20, height = 14, units = "cm")

## 3. Spatial model ----
names(van.island.UTM)
van.island.UTM

# split out UTM X and Y columns
utm.coords <- data.frame(st_coordinates(van.island.UTM))
utm.coords <- utm.coords %>% 
  mutate(X = X / 1000,
         Y = Y / 1000) # convert from metres to kilometres
head(utm.coords)
van.island.UTM <- cbind(van.island.UTM, utm.coords)

# define Loma salmonae presence/absence column

van.island.UTM$lo_sal_binary <- 0
van.island.UTM$lo_sal_binary[van.island.UTM$lo_sal > 0] <- 1
table(van.island.UTM$lo_sal_binary) # Loma salmonae found in 1229 of 5060 samples

ggplot() +
  geom_sf(data = van.island.UTM, aes(colour = factor(lo_sal_binary))) +
  facet_wrap(~ Year)

# make an sdmTMB model

dat <- data.frame(van.island.UTM)
class(dat) # sdmTMB requires a data.frame, not an sf object

loma.mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)
plot(loma.mesh)

# fit a very simple model: 
# Day of year: hypothesize Loma salmonae infections to increase with more time
# in the ocean, since this is a marine-transmitted parasite
dat$doy <- yday(dat$Date)
hist(dat$doy)

# Population: include a random intercept to account for inter-population variation
dat$Stock_Region <- factor(dat$Stock_Region)

m1 <- sdmTMB(lo_sal_binary ~ s(doy, bs = "cc") + (1 | Stock_Region),
             data = dat,
             mesh = loma.mesh,
             family = binomial(link = "logit"), # binomial GLM since presence/absence
             spatial = "on")
m1
sanity(m1) # provides a method to check model

# before we look at what the model says, how well does it fit? Let's look at the
# residuals

# QQ plot
m1.res <- residuals(m1, type = "mle-mvn")
qqnorm(m1.res);abline(0, 1) # good

# histogram of residuals
hist(m1.res) # good

# plot residuals over space: no obvious spatial pattern, which is good
dat$residual <- residuals(m1, type = "mle-mvn")

ggplot() +
  geom_point(data = dat,
             aes(x = X, y = Y, colour = residual), alpha = 0.6) +
  scale_color_distiller(palette = "RdBu") +
  facet_wrap(~ Year) 

# plot residuals against covariates
ggplot(data = dat,
       aes(x = doy, y = residual)) +
  geom_point(alpha = 0.3) +
  geom_smooth() #no obvious temporal pattern - good

ggplot(dat = dat,
       aes(x = Stock_Region, y = residual)) +
  geom_jitter(alpha = 0.3) +
  geom_smooth() +
  theme(axis.text.x = element_text(angle = 90)) # no obvious pattern among populations - good
# Conclusion: model fits pretty well!

## to visualize the model, we need to make a prediction grid!
# We will use "conditional predictions," where the model predicts across the 
# range of the variable of interest, while other variables are fixed at their 
# median (quantitative variables) or most common (categorical variables) value.

grid.4km <- st_make_grid(van.island.UTM, cellsize = 4000, # 2km grid
  crs = st_crs(32610), what = "centers")
grid.4km <- st_as_sf(grid.4km)

# identify which grid cells are on land and remove them
grid.4km.join <- st_join(x = grid.4km, y = bc.coast.UTM)
head(grid.4km.join)

grid.4km.ocean <- subset(grid.4km.join, is.na(NAME))
class(grid.4km.ocean)

ggplot() +
  geom_sf(data = grid.4km.ocean) +
  geom_sf(data = van.island.UTM, colour = "red")
# there are so many points that are very far from sampling areas,
# so let's remove those

# calculate distance to shore
sw.bc <- readRDS("C:/Users/wgree/OneDrive/Documents/ASDP_biogeog/_Nathanael/ASDP_Index/data/sw-bc.RDS")
st_crs(sw.bc)
sw.bc <- st_transform(sw.bc, crs = 32610)

shore.distance <- st_distance(x = grid.4km.ocean, y = sw.bc)
grid.4km.ocean$shore_dist <- as.numeric(shore.distance[, 1])

ggplot() +
  geom_sf(data = subset(grid.4km.ocean, shore_dist < 60000), 
          aes(colour = shore_dist)) +
  geom_sf(data = van.island.UTM, colour = "red") +
  scale_colour_viridis_c(option = "F")

grid.4km.ocean <- subset(grid.4km.ocean, shore_dist < 60000)
grid.coords <- data.frame(st_coordinates(grid.4km.ocean))

spatial.predict.grid <- expand.grid(
  coord = paste(grid.coords$X / 1000,
                grid.coords$Y / 1000, sep = ","), # spatial grid
  doy = 183,
  Stock_Region = factor("ECVI"))

# divide spatial.predict.grid$coord into X and Y
spatial.predict.grid <- separate(spatial.predict.grid, coord, c("X", "Y"), 
                                 sep = ",", remove = TRUE)
spatial.predict.grid$X <- as.numeric(spatial.predict.grid$X)
spatial.predict.grid$Y <- as.numeric(spatial.predict.grid$Y)

m1.spatial.pred <- predict(m1, newdata = spatial.predict.grid, 
                           type = "response", se_fit = FALSE)
hist(m1.spatial.pred$est)
ggplot() +
  geom_raster(data = m1.spatial.pred, 
            aes(X * 1000, Y * 1000, fill = est)) +
  geom_sf(data = bc.coast.UTM) +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_sf(xlim = c(70 * 1000, 525 * 1000), 
           ylim = c(5320 * 1000, 5680 * 1000)) +
  scale_x_continuous(breaks = c(-129, -127, -125, -123)) +
  scale_y_continuous(breaks = c(48, 49, 50, 51)) +
  theme(axis.title = element_blank()) +
  labs(fill = "Prevalence")
