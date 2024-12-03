#### script to identify rivers and create nearest distance to water raster

#### set up environemnt ####
install.packages("whitebox", repos="http://R-Forge.R-project.org")
library(whitebox)
library(terra) 
library(rsi)
library(sf)
library(fs)
library(raster)
library(stars)

####use water shapefiles ####


# Define file paths
rivers_raster <- rast("C:/workspace/PhD year 2/Thesis_work/Data/HydroRiver_raster.tif")
reference_raster <- rast("Data/combined_metrics_raster.tif")
output_distance_raster_path <- "Data/Outputs/distance_to_water.tif"

rivers_raster <- classify(rivers_raster, cbind(-Inf, 0, NA))
rivers_raster[rivers_raster > 0] <- 1   
# Calculate Euclidean distance from each pixel to the nearest road
distance_to_rivers <- distance(rivers_raster, target = NA)

# Clip the distance raster to match the extent of the reference raster
distance_to_rivers_clipped <- crop(distance_to_rivers, reference_raster)

# Save the result to a new file
writeRaster(distance_to_rivers_clipped, output_distance_raster_path, overwrite = TRUE)
print(paste("Distance to rivers raster saved to:", output_distance_raster_path))

