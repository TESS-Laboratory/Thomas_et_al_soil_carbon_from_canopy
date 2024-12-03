#### script to stack layers and extract data for models

### set up envoironment ####
library(terra)
library(raster)
library(rsi)

#### read in data ####

##add in function to read in lidar raster stack

folder_path <-"Data/SoilGrids"  # Replace with the actual folder path
SG_files <- list.files(path = folder_path, pattern = "\\.tif$", full.names = TRUE)

SG_first<- rast(SG_files[1])
SG_list <- lapply(SG_files, rast)
SG_aligned <- lapply(SG_list, resample, y= SG_first, method = "bilinear")
SG<- rast(SG_aligned)

output<- "Data/output"

#### merge rasters ####
raster_list <- lapply(stack_list, rast)  # Read all rasters as stacks
M<- do.call(mosaic, raster_list)


mosaic_layers<- rast(M)

# Align CRS
crs(SG) <- crs(mosaic_layers)

# Resample target raster to match resolution and extent of reference raster
aligned_raster <- resample(mosaic_layers, SG, method = "bilinear")

matrics<- c(SG, aligned_raster)
#metrics<- rsi::stack_rasters(list(SG, mosaic_layers), output_filename = paste0(output, "/metricstack.vrt"))


#### debugging ####
raster<-raster_list[[3]]
layer <- raster[[2]]
band_name <- paste0(names(layer))

# If this is the first time encountering this band name, create an entry in the list
if (!band_name %in% names(band_layers)) {
  band_layers[[band_name]] <- layer
} else {
    band_layers[[band_name]] <- mosaic(band_layers[[band_name]], layer[[band_name]])
}


R<- mosaic(raster_list[[1]], raster_list[[2]])
