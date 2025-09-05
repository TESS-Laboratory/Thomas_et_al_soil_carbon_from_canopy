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

### set up extraction parameters ####
fp<- "C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data"
metrics<- rast(file.apth(fp, "combined_metrics_raster.tif"))
samples<- read_csv(file,path(fp, "soil_meta_table.csv"))
coordinates <- vect(samples[, c("plot.x", "plot.y", "Codigo")], geom = c("plot.x", "plot.y"), crs = "EPSG:4326")

#align coordinates to rasters
proj<- crs(metrics)
coordinates <- project(coordinates, proj, partial = TRUE)

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

gedi_2B_sf <- grab_gedi(gedi_2B_search) |>
  filter(
    l2a_quality_flag == 1,
    l2b_quality_flag == 1,
    degrade_flag == 0,
    dplyr::between(sensitivity, 0.98, 1.0) # upped from 0.9 to 0.98
  ) |>
  dplyr::select(
    shot_number,
    date_time,
    cover,
    fhd_normal,
    pai,
    pgap_theta,
    rv,
    dplyr::starts_with("cover_z"),
    dplyr::starts_with("pai_z"),
    dplyr::starts_with("pavd_z"),
    lat_lowestmode,
    lon_lowestmode,
  ) |>
  collect_gedi(gedi_find = gedi_2B_search)


chewie_show(
  gedi_2B_sf,
  zcol = "pai",
  zoom = 8,
  alpha = 0.5,
  aoi_color = "white"
)
Gedi_repro<- st_transform(gedi_2B_sf, crs = proj)
extracted_values <- st_join(coordinates, Gedi_repro,join= st_nearest_feature, bind = TRUE)
st_write(extracted_values, file.path(fp, "Gedi_2b_dataframe.csv"))
