# Interative maps

library(leaflet)
library(sf)

rockfish <- read.csv("data/iNaturalist-rockfish.csv")

# make a map to start
## |> is the baseR equivalent of the pipe ( %>% )

leaflet() |>
  setView(lng = -123, lat = 48.4, zoom = 9) |>
  addProviderTiles(providers$Esri.WorldImagery) |>
  addCircles(data = rockfish, 
             lng = rockfish$longitude, lat = rockfish$latitude,
             color = "red")
