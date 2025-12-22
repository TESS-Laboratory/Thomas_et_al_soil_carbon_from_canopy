#### script to stack layers and extract data for models
## make sure all spatial data is in the same resolution or add terra::project to make it so before merging
#this was performed with epsg:31980
### set up envoironment ####
install.packages(c("rsi", "terra", "patchwork", "ggplot2", "viridis"), repos = "https://cloud.r-project.org")
library(terra)
library(rsi)
library(patchwork)
library(ggplot2)
library(viridis)

fp<-"~/workspace/PhD_work/soil_chapter/Data" ##workstation
#fp<- "C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data"
#### read in data ####
vrt.folder <- file.path("~/workspace/PhD_work/Thomas_et_al_soil_carbon_from_canopy/Data/RC") ##replace with location of vrt files computed by lidar_analysis.R
tif.folder <- file.path(fp, "SoilGrids") ## Replace with the location of Tiff files (soil grids, rivers, ect)
output<- file.path("~/workspace/PhD_work/Thomas_et_al_soil_carbon_from_canopy/Data") ## replace with desired location of output files


##add in function to read in lidar raster stack
# Function to list all files in each subfolder of a specified folder
list_files_by_subfolder <- function(folder_path) {
  # Get immediate subfolders
  subfolders <- list.dirs(folder_path, full.names = TRUE, recursive = TRUE)
  
  # Initialize output list
  files_by_subfolder <- list()
  
  for (subfolder in subfolders) {
    
    # Find .tif files in this subfolder
    files <- list.files(
      subfolder,
      full.names = TRUE,
      pattern = "\\.tif$"
    )
    
    # Skip folders with no matching files
    if (length(files) == 0) next
    
    # Use folder name (not undefined 'folder') as list name
    rasters<- list()
    for (file in files) {
      R<- rast(file)
      rasters[[file]]<- R
    }
    
    folder_name <- basename(subfolder)
    
    # Read files as SpatRaster
    files_by_subfolder[[folder_name]] <- mosaic(sprc(rasters))
  }
  
  return(files_by_subfolder)
}


# Replace with the actual folder path
nested_files_list <- list_files_by_subfolder(vrt.folder)
#names(nested_files_list)<- NULL
#### merge rasters ####
res<- nested_files_list[[1]]
lidar_aligned <- lapply(nested_files_list, terra::resample, y= res, method = "bilinear")
raster_list <- sprc(lidar_aligned)  # Read all rasters as stacks


mos<- mosaic(raster_list, overwrite = TRUE)



SG_files <- list.files(path = tif.folder, pattern = "\\.tif$", full.names = TRUE)

SG_first<- rast(SG_files[2])
SG_list <- lapply(SG_files, rast)
names(SG_list[[11]])<- paste0("year_of_last_fire")
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

writeRaster( metrics, file.path(fp, "combined_metrics_raster.tif"), overwrite = TRUE)
