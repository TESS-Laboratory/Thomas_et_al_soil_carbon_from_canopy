### script to obtain GEDI data via chewie
### set up environment ###
# install.packages("pak")
pak::pkg_install("Permian-Global-Research/chewie")
install.packages('mapview')
library(chewie)
library(dplyr)
library(sf)
library(mapview)
chewie_creds() # to set up your credentials
chewie_health_check() # to check your credentials and cache setup.

## search for gedi data ####
RC_Boundary <- sf::st_read("C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Rio_Cautario/RC_boundary_EPSG4326/RC_boundary_EPSG4326.shp", crs = 4326
)

gedi_2B_search <- find_gedi(RC_Boundary,
                            gedi_product = "2B",
                            date_start = "2022-05-01",
                            date_end = "2024-07-31"
)
print(gedi_2B_search)

chewie_show(
  gedi_2B_search,
  zoom = 8
)

gedi_2B_sf <- grab_gedi(gedi_2a_search) |>
  filter(
    quality_flag == 1,
    degrade_flag == 0
  ) |>
  select(
    beam, date_time, lat_lowestmode, lon_lowestmode, elev_highestreturn,
    elev_lowestmode, rh0, rh25, rh50, rh75, rh95, rh100
  ) |>
  collect_gedi(gedi_find = gedi_2a_search)
