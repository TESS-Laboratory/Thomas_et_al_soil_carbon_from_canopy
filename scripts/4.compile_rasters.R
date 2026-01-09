#### script to stack layers and extract data for models
## make sure all spatial data is in the same resolution or add terra::project to make it so before merging
#this was performed with epsg:31980
### set up envoironment ####
#install.packages("rsi")
library(terra)
library(rsi)
library(patchwork)
library(ggplot2)
library(viridis)


fp<- "C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data"
#### read in data ####

tif.folder <- file.path(fp, "SoilGrids") ## Replace with the location of Tiff files (soil grids, rivers, ect)
mosaic_layers<- rast(file.path(fp, "lidar_metrics.tif"))


SG_files <- list.files(path = tif.folder, pattern = "\\.tif$", full.names = TRUE)

SG_first<- rast(SG_files[2])
SG_list <- lapply(SG_files, rast)
names(SG_list[[4]])<- paste0("min_distance")
names(SG_list[[17]])<- paste0("year_of_last_fire")
SG_aligned <- lapply(SG_list, terra::resample, y= SG_first, method = "bilinear")
SG<- rast(SG_aligned)
names(SG) <- paste0(names(SG), "_1")




#project(SG_first, "epsg:31980")
# Align CRS
crs(SG) <- crs(mosaic_layers)
names(mosaic_layers) <- paste0(names(mosaic_layers), "_3")

# Resample target raster to match resolution and extent of reference raster
aligned_raster <- terra::resample(mosaic_layers, SG, method = "bilinear")

metrics<- c(SG, aligned_raster)

writeRaster( metrics, file.path(fp, "combined_metrics_raster.tif"), overwrite = TRUE)
