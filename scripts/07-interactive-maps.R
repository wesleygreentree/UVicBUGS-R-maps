# Interative maps
# Note |> is the baseR equivalent of the pipe ( %>% ) used in the tidyverse

library(leaflet)
library(sf)
library(tidyverse)

rockfish <- read.csv("data/iNaturalist-rockfish.csv")

# standardize capitalization of species name
rockfish$common_name <- tolower(rockfish$common_name)

# make a legend for species
length(unique(rockfish$common_name)) # there are 25 rockfish species in the dataset,
# which is too many for a legend. Simplify:

species.tally <- data.frame(table(rockfish$common_name))
common.species <- species.tally |> arrange(desc(Freq)) |> head(n = 7)

rockfish$species_legend <- rockfish$common_name
rockfish$species_legend[!rockfish$common_name %in% common.species$Var1] <- "other"
table(rockfish$species_legend)

rockfish$species_legend <- factor(rockfish$species_legend, 
  levels = c("black rockfish", "brown rockfish", "copper rockfish", "quillback rockfish",
             "tiger rockfish", "yelloweye rockfish", "yellowtail rockfish", "other"))

# make colour palette
rockfish.pal <- colorFactor(palette = "Set2", levels(rockfish$species_legend), ordered = TRUE)

# make label for each observation
labels <- sprintf(
  "<strong>Species: </strong>%s<br/><strong>Date: </strong>%s",
  rockfish$common_name, rockfish$observed_on) %>% 
  lapply(htmltools::HTML)

leaflet() |>
  setView(lng = -124, lat = 48.8, zoom = 7) |>
  addProviderTiles(providers$Esri.WorldImagery) |>
  addCircles(data = rockfish, 
             lng = rockfish$longitude, lat = rockfish$latitude,
             color = ~rockfish.pal(rockfish$species_legend), 
             label = labels, opacity = 1) |>
  addLegend(pal = rockfish.pal, values = levels(rockfish$species_legend), opacity = 0.7, 
            title = "Rockfish species", position = "topright")

#Helpful resource: https://leafletjs.com/examples/choropleth/