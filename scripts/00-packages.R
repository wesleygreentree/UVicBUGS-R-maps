# Install required packages
# if you have some but not all, feel free to only install those you don't
# have currently

install.packages(c("sf", "stars", "tidyverse",
                   "basemaps", "rphylopic", "gganimate",
                   "terra", "tidyterra", "palmerpenguins",
                   "rnaturalearth", "rnaturalearthdata", "leaflet", "ggnewscale",
                   "patchwork", "mapedit", "ggspatial"))

install.packages("rnaturalearthhires",
  repos = "https://ropensci.r-universe.dev",
  type = "source") # high resolution basemaps for RNaturalEarth

