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

## 2C. Calculate prevalence in each grid cell  ----

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

## 2D Leaflet interactive map ----
lo.sal.LL <- st_transform(lo.sal, crs = st_crs(4326)) # return shapefile to WGS 84

# make colour scale
bins <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
pal <- colorBin("YlOrRd", domain = lo.sal.LL$Percentage, bins = bins)
# uses RColorBrewer palettes: colorbrewer2.org

# make labels
labels <- sprintf(
  "<strong>Total sample size: </strong>%s<br/>L. salmonae: %g",
  lo.sal.LL$total_n, round(lo.sal.LL$Percentage, 1)) %>% 
  lapply(htmltools::HTML)

leaflet(lo.sal.LL) %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>% 
  addPolygons(fillColor = ~pal(Percentage),
              label = labels,
              fillOpacity = 0.8, stroke = FALSE) %>% 
  addLegend(pal = pal, values = ~Percentage, opacity = 0.7, title = "Percentage",
            position = "topright")
  
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

m1 <- sdmTMB(lo_sal_binary ~ s(doy, bs = "cc", k = 5) + (1 | Stock_Region),
             data = dat,
             mesh = loma.mesh,
             family = binomial(link = "logit"), # binomial GLM since presence/absence
             spatial = "on",
             knots = list(doy = c(1, 365))) # first and last day of year, needed
             # for the cyclic cubic spline
m1
sanity(m1) # provides a method to check model
AIC(m1) # 4779.6

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
sw.bc <- readRDS("data/sw-bc.RDS")
st_crs(sw.bc)
sw.bc <- st_transform(sw.bc, crs = 32610)

shore.distance <- st_distance(x = grid.4km.ocean, y = sw.bc) # takes a while
grid.4km.ocean$shore_dist <- as.numeric(shore.distance[, 1])

ggplot() +
  geom_sf(data = subset(grid.4km.ocean, shore_dist < 60000), 
          aes(colour = shore_dist)) +
  geom_sf(data = van.island.UTM, colour = "red") +
  scale_colour_viridis_c(option = "F")

grid.4km.ocean <- subset(grid.4km.ocean, shore_dist < 60000)
grid.coords <- data.frame(st_coordinates(grid.4km.ocean))

# predict spatially
spatial.predict.grid <- expand.grid(
  coord = paste(grid.coords$X / 1000,
                grid.coords$Y / 1000, sep = ","), # spatial grid
  doy = 183, # median(dat$doy) = 183
  Stock_Region = factor("ECVI"))

# divide spatial.predict.grid$coord into X and Y
spatial.predict.grid <- separate(spatial.predict.grid, coord, c("X", "Y"), 
                                 sep = ",", remove = TRUE)
spatial.predict.grid$X <- as.numeric(spatial.predict.grid$X)
spatial.predict.grid$Y <- as.numeric(spatial.predict.grid$Y)

m1.spatial.pred <- predict(m1, newdata = spatial.predict.grid, 
                           type = "response", se_fit = FALSE)

ggplot() +
  geom_raster(data = m1.spatial.pred, 
            aes(X * 1000, Y * 1000, fill = est)) +
  geom_sf(data = bc.coast.UTM) +
  scale_fill_viridis_c(option = "B", trans = "identity") +
  coord_sf(xlim = c(70 * 1000, 525 * 1000), 
           ylim = c(5320 * 1000, 5680 * 1000)) +
  scale_x_continuous(breaks = c(-129, -127, -125, -123)) +
  scale_y_continuous(breaks = c(48, 49, 50, 51)) +
  theme(axis.title = element_blank()) +
  labs(fill = "Estimated\nprevalence")

m1.spatial.pred.uncertainty <- predict(m1, newdata = spatial.predict.grid, 
                           type = "response", se_fit = TRUE)
# returns predictions on link-scale
# the predictions are logit-transformed, we exponentiate to bring them back
# to the original scale
m1.spatial.pred.uncertainty$se_response <- exp(m1.spatial.pred.uncertainty$est_se)
hist(m1.spatial.pred.uncertainty$se_response)

ggplot() +
  geom_raster(data = m1.spatial.pred.uncertainty, 
              aes(X * 1000, Y * 1000, fill = se_response)) +
  geom_sf(data = bc.coast.UTM) +
  scale_fill_viridis_c(trans = "identity") +
  coord_sf(xlim = c(70 * 1000, 525 * 1000), 
           ylim = c(5320 * 1000, 5680 * 1000)) +
  scale_x_continuous(breaks = c(-129, -127, -125, -123)) +
  scale_y_continuous(breaks = c(48, 49, 50, 51)) +
  theme(axis.title = element_blank()) +
  labs(fill = "Standard\nerror")

## spatial predictions for each season
ggplot() +
  geom_jitter(data = dat,
             aes(x = doy, y = 1, colour = season2),
             width = 0, height = 1) +
  xlim(0, 365)

# we will take the median sample day per season, but there are many ways you could do this!
dat %>% 
  mutate(doy = yday(Date)) %>% 
  group_by(season2) %>% 
  summarize(median_doy = median(doy))

seasonal.predict.grid <- expand.grid(
  coord = paste(grid.coords$X / 1000,
                grid.coords$Y / 1000, sep = ","), # spatial grid
  doy = c(169, 264),
  Stock_Region = factor("ECVI"))

# divide seasonal.predict.grid$coord into X and Y
seasonal.predict.grid <- separate(seasonal.predict.grid, coord, c("X", "Y"), 
                                 sep = ",", remove = TRUE)
seasonal.predict.grid$X <- as.numeric(seasonal.predict.grid$X)
seasonal.predict.grid$Y <- as.numeric(seasonal.predict.grid$Y)

m1.seasonal.pred <- predict(m1, newdata = seasonal.predict.grid, 
                           type = "response", se_fit = FALSE)

m1.seasonal.pred <- m1.seasonal.pred %>% 
  mutate(label = recode(as.character(doy),
                        "169" = "June", "264" = "September"))

ggplot() +
  geom_raster(data = m1.seasonal.pred, 
              aes(X * 1000, Y * 1000, fill = est)) +
  geom_sf(data = bc.coast.UTM) +
  scale_fill_viridis_c(option = "B", trans = "identity") +
  coord_sf(xlim = c(70 * 1000, 525 * 1000), 
           ylim = c(5320 * 1000, 5680 * 1000)) +
  scale_x_continuous(breaks = c(-129, -127, -125, -123)) +
  scale_y_continuous(breaks = c(48, 49, 50, 51)) +
  theme(axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11)) +
  labs(fill = "Estimated\nprevalence") +
  facet_wrap(~ label)

# Last, plot the fit of the GAM smooth from svc1.seasonal.pred

# to account for spatial variation, let's predict the smooth at the four
# most-sampled grids

most.common.grids <- grid.sample.size %>% arrange(desc(n)) %>% head(n = 4)

ggplot() +
  geom_sf(data = bc.coast.UTM) +
  geom_sf(data = subset(hexagon.grid, grid_id %in% most.common.grids$grid_id),
          fill = "red") +
  coord_sf(xlim = c(70 * 1000, 525 * 1000), 
           ylim = c(5320 * 1000, 5680 * 1000))

# take the grid centroids
most.common.grids.centroids <- data.frame(st_coordinates(st_centroid(most.common.grids)))
names(most.common.grids.centroids)

doy.predict.grid <- expand.grid(
  coord = paste(most.common.grids.centroids$X / 1000, 
                most.common.grids.centroids$Y / 1000, sep = ","),
  doy = seq(1, 365, by = 1), # every third day
  Stock_Region = factor("ECVI") # most common stock
)

# divide doy.grid$coord into X and Y
doy.predict.grid <- separate(doy.predict.grid, coord, c("X", "Y"), 
                             sep = ",", remove = TRUE)
doy.predict.grid$X <- as.numeric(doy.predict.grid$X)
doy.predict.grid$Y <- as.numeric(doy.predict.grid$Y)

# assign season2: April to August, September to March
table(dat$season2)
doy.predict.grid$season2 <- "FaWi"
doy.predict.grid$season2[doy.predict.grid$doy %in% 91:243] <- "SpSu"

m1.doy.pred <- predict(m1, newdata = doy.predict.grid, 
                         type = "response", se_fit = TRUE)
# returns predictions at link scale

# get confidence intervals
m1.doy.pred$lo_ci <- m1.doy.pred$est + (qnorm(0.025) * m1.doy.pred$est_se)
m1.doy.pred$up_ci <- m1.doy.pred$est + (qnorm(0.975) * m1.doy.pred$est_se)

# return to response scale with the inverse logit: exp(x)/(1+exp(x))
m1.doy.pred$est_response <- exp(m1.doy.pred$est) / (1 + exp(m1.doy.pred$est))
m1.doy.pred$lo_ci_response <- exp(m1.doy.pred$lo_ci) / (1 + exp(m1.doy.pred$lo_ci))
m1.doy.pred$up_ci_response <- exp(m1.doy.pred$up_ci) / (1 + exp(m1.doy.pred$up_ci))

# plot one smooth per coordinate
m1.doy.pred$coord <- paste(m1.doy.pred$X, m1.doy.pred$Y, sep = "-")

ggplot() +
  geom_ribbon(data = m1.doy.pred,
              aes(x = doy, ymin = lo_ci_response, ymax = up_ci_response,
                  fill = coord), 
              alpha = 0.3, show.legend = FALSE) +
  geom_line(data = m1.doy.pred,
            aes(x = doy, y = est_response, colour = coord),
            show.legend = FALSE)
# which site does each smooth represent
ggplot() +
  geom_sf(data = bc.coast.UTM) +
  geom_point(data = m1.doy.pred,
             aes(x = X * 1000, y = Y * 1000, colour = coord),
             size = 3, show.legend = FALSE) +
  coord_sf(xlim = c(70 * 1000, 525 * 1000), 
           ylim = c(5320 * 1000, 5680 * 1000))
            
## 3B Spatially-varying coefficients ----
## take a look at spatially-varying coefficients, which is where the latent spatial
# variation can vary by a covariate
# Here we will let the spatial variation vary by season
dat$season2 <- factor(dat$season2)
unique(dat$season2) # Spring summer, fall winter
class(dat$season2)

svc1 <- sdmTMB(lo_sal_binary ~ s(doy, bs = "cc", k = 5) + (1 | Stock_Region),
               data = dat,
               mesh = loma.mesh,
               spatial_varying = ~ 0 + season2,
               family = binomial(link = "logit"), # binomial GLM since presence/absence
               spatial = "off",
               knots = list(doy = c(1, 365))) # spatial_varying ~ 0 + season2 includes the spatial
                                # intercept so spatial = "on" is not needed
#Use the Akaike information criterion to evaluate model fit
AIC(m1)
AIC(svc1) # lower AIC value indicates that svc1 is better than m1
sanity(svc1)

# quick look at residuals
svc1.res <- residuals(svc1, type = "mle-mvn")
qqnorm(svc1.res);abline(0, 1) # good

hist(svc1.res) # good

# plot residuals over space:
dat$svc_residual <- residuals(svc1, type = "mle-mvn")

ggplot() +
  geom_point(data = dat,
             aes(x = X, y = Y, colour = svc_residual), alpha = 0.6) +
  scale_color_distiller(palette = "RdBu") +
  facet_wrap(~ Year) 
# may be worth adding a covariate to see if it improves residauls in Strait of 
# Georgia in a couple years, but fairly acceptable

# predict from model
seasonal.predict.grid <- seasonal.predict.grid %>% 
  mutate(season2 = recode(as.character(doy),
                         "169" = "SpSu", "264" = "FaWi"))

svc1.seasonal.pred <- predict(svc1, newdata = seasonal.predict.grid, 
                            type = "response", se_fit = FALSE)

svc1.seasonal.pred <- svc1.seasonal.pred %>% 
  mutate(label = recode(season2,
                        "SpSu" = "June", "FaWi" = "September"))

ggplot() +
  geom_raster(data = svc1.seasonal.pred, 
              aes(X * 1000, Y * 1000, fill = est)) +
  geom_sf(data = bc.coast.UTM) +
  scale_fill_viridis_c(option = "B", trans = "identity") +
  coord_sf(xlim = c(70 * 1000, 525 * 1000), 
           ylim = c(5320 * 1000, 5680 * 1000)) +
  scale_x_continuous(breaks = c(-129, -127, -125, -123)) +
  scale_y_continuous(breaks = c(48, 49, 50, 51)) +
  theme(axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11)) +
  labs(fill = "Estimated\nprevalence") +
  facet_wrap(~ label)
