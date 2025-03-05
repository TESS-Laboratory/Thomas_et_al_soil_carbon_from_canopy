#### script to stack layers and extract data for models
## make sure all spatial data is in the same resolution or add terra::project to make it so before merging
#this was performed with epsg:31980
### set up envoironment ####
library(terra)
library(rsi)



#### read in data ####
vrt.folder <- "C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Large_Data/RC" ##replace with location of vrt files computed by lidar_analysis.R
tif.folder <-"C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data/SoilGrids"  ## Replace with the location of Tiff files (soil grids, rivers, ect)
output<- "Data/Outputs" ## replace with desired location of output files


##add in function to read in lidar raster stack
# Function to list all files in each subfolder of a specified folder
list_files_by_subfolder <- function(folder_path) {
  # Get a list of all subfolders within the main folder
  subfolders <- list.dirs(folder_path, full.names = TRUE, recursive = FALSE)
  
  # Initialize a list to store file paths by subfolder
  files_by_subfolder <- list()
  
  # Iterate over each subfolder
  for (subfolder in subfolders) {
    # Get a list of all files in the current subfolder
    files <- list.files(subfolder, full.names = TRUE, pattern = "\\comp.vrt")
    
    # Store the files in the list, named by the subfolder path
    files_by_subfolder[[subfolder]] <- rast(files)
  }
  
  # Return the list of files grouped by subfolder
  return(files_by_subfolder)
}


# Replace with the actual folder path
nested_files_list <- list_files_by_subfolder(vrt.folder)
names(nested_files_list)<- NULL
#### merge rasters ####
raster_list <- sprc(nested_files_list)  # Read all rasters as stacks
mos<- mosaic(raster_list, overwrite = TRUE)



SG_files <- list.files(path = tif.folder, pattern = "\\.tif$", full.names = TRUE)

SG_first<- rast(SG_files[2])
SG_list <- lapply(SG_files, rast)
SG_aligned <- lapply(SG_list, terra::resample, y= SG_first, method = "bilinear")
SG<- rast(SG_aligned)
names(SG) <- paste0(names(SG), "_1")




mosaic_layers<- mos
#project(SG_first, "epsg:31980")
# Align CRS
crs(SG) <- crs(mosaic_layers)
names(mosaic_layers) <- paste0(names(mosaic_layers), "_3")

# Resample target raster to match resolution and extent of reference raster
aligned_raster <- terra::resample(mosaic_layers, SG, method = "bilinear")

metrics<- c(SG, aligned_raster)

writeRaster( metrics, "C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Large_Data/combined_metrics_raster.tif", overwrite = TRUE)
#### debugging ####
