
library(tidyverse)
library(sf)
library(readxl)

bc.coast <- read_sf("data/bc-coast.shp")

chinook.disease <- read_excel("data/disease-data.xlsx",
                              guess_max = 5719, sheet = "chinook")
ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = chinook.disease,
             aes(x = Longitude, y = Latitude), alpha = 0.3) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58)) +
  theme_bw()

table(chinook.disease$Stock_Region)

ggplot() +
  geom_sf(data = bc.coast) +
  geom_point(data = subset(chinook.disease, Stock_Region == "ECVI"),
             aes(x = Longitude, y = Latitude), alpha = 0.3) +
  coord_sf(xlim = c(-137, -122), ylim = c(47, 58)) +
  theme_bw()
view(subset(chinook.disease, Stock_Region == "ECVI"))
