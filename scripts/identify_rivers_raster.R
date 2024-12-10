#### script to identify rivers and create nearest distance to water raster

#### set up environemnt ####
library(terra) 
#library(rsi)
#library(sf)
#library(fs)


####use water shapefiles ####


# Define file paths
rivers_raster <- rast("Data/HydroRiver_raster.tif") ## file path to binarised river tif
reference_raster <- rast("Data/combined_metrics_raster.tif") ## file path to raster in CRS and resolution you need
final_tif_ouput<-"Data/SoilGrids/distance_to_water.tif"## file path you want the output to be saved (ideally with other tifs that have been downloaded)

rivers_raster <- classify(rivers_raster, cbind(-Inf, 0, NA))
rivers_raster[rivers_raster > 0] <- 1   
# Calculate Euclidean distance from each pixel to the nearest road
distance_to_rivers <- distance(rivers_raster, target = NA)

# Clip the distance raster to match the extent of the reference raster
distance_to_rivers_clipped <- crop(distance_to_rivers, reference_raster)

# Save the result to a new file
writeRaster(distance_to_rivers_clipped, final_tif_ouput, overwrite = TRUE)
print(paste("Distance to rivers raster saved to:", final_tif_ouput))

