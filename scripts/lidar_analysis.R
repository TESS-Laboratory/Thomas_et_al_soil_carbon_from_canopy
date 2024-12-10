####Lidar analysis script to get canopy structure metrics from aiborne lidar data  
#### best run on workstation ####


#### set up environment ####
install.packages("lidR")
install.packages("filesstrings")
install.packages("SpaDES")
install.packages("terra")
install.packages("geometry")
install.packages("purrr")
install.packages("rsi")

library(lidR)
library(geometry)
#library(raster)
library(terra)
library(stars)
library(sp)
library(sf)
library(tidyverse)
library(rlas)
library(rgl)
library(filesstrings)
library(SpaDES)
library(purrr)
library(rsi)

future::plan(future::multisession, workers = future::availableCores()-2)

###### read in functions ######

#' @title Create a file if it does or does not exist.
#' @description Creates a title and prints useful information
create_my_dir <- function(dir, purpose = "") {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cli::cli_alert_info(
      "Created {purpose}  Directory in : '{dir}'"
    )
  } else {
    cli::cli_alert_info(
      "{purpose} Directory already exists in : '{dir}'"
    )
  }
  return(dir)
}

#' @title setup {lidR} LAScatalog object
#' @description setup {lidR} LAScatalog object with opinionated vars and
#' file paths
#' @param laz filepath for laz file.
#' @param area_name name of area
#' @param proj_name name of project area
#' @param crs projection of the lidar data, numeric epsg code or character
#' proj4string is fine.
#' @param .force force the laz file to be copied to the processing directory
#' @return a LAScatalog object
lidar_setup_ctg <- function(
    laz, area_name, proj_name,
    chunk_size = 260,
    chunk_buffer = 10,
    crs = NULL,
    .force = FALSE,
    proc_dir_root = here::here("Data")) {
 if(!dir.exists(proc_dir_root)) dir.create(proc_dir_root)
   lidar_processing_dir <- create_my_dir(
    file.path(proc_dir_root, proj_name, area_name)
  )
  cp_dest <- file.path(
    lidar_processing_dir, basename(laz)
  )
  
  if (!file.exists(cp_dest) || .force) {
    cli::cli_alert_info("Copying LiDAR data to processing directory")
    file.copy(laz, cp_dest)
  } else {
    cli::cli_alert_info("LiDAR data already exists in processing directory")
  }
  
  laz_ctg <- lidR::readLAScatalog(dirname(cp_dest))
  
  # set options
  if (!is.null(chunk_size)) {
    lidR::opt_chunk_size(laz_ctg) <- chunk_size
  }
  if (!is.null(chunk_buffer)) {
    lidR::opt_chunk_buffer(laz_ctg) <- chunk_buffer
  }
  if (!is.null(crs)) {
    lidR::st_crs(laz_ctg) <- crs
  }
  
  lidR::opt_laz_compression(laz_ctg) <- TRUE
  
  
  return(laz_ctg)
}


## list all lidar files
las_dir<-"C:/workspace/PhD year 2/datasets/Rio_Cautario_ALS/LiDAR"
file.name<- grep(list.files(path=las_dir, full.names = TRUE), pattern= '*.copc.laz', invert=TRUE, value=TRUE)

##set up and run function across all files
area_names<- basename(file.name)%>% 
  tools::file_path_sans_ext()

catalog_list<- purrr::map2(
  .x= file.name, .y=area_names, 
  ~lidar_setup_ctg(.x, .y, "RC", crs= "EPSG:31980" ))

##### normalise point cloud and create dtm ####

#' @title Generate a DTM from a laz catalog
#' @description Generate a DTM from a laz catalog with opinoated settings
#' @param laz_ctg A LAScatalog object
#' @param .alg A function to use for generating the DTM
#' @param .force Force the DTM to be generated even if it already exists
#' @return A SpatRaster object
lidar_build_dtm <- function(
    laz_ctg,
    .alg = lidR::knnidw(5L),
    .force = FALSE,
    .return_class = c("SpatRaster", "character")) {
  checkmate::assert_class(laz_ctg, "LAScatalog")
  checkmate::assert_class(.alg, "function")
  checkmate::assert_class(.force, "logical")
  .return_class <- rlang::arg_match(.return_class)
  lid_proc_dir <- dirname(laz_ctg@data$filename[1])
  mod_name <- basename(lid_proc_dir)
  proc_dir <- file.path(lid_proc_dir, paste0(mod_name, "_dtm_tiles"))
  lidR::opt_output_files(laz_ctg) <- paste(
    proc_dir,
    paste0(mod_name, "_dtm_{ID}"),
    sep = "/"
  )
  vrt_check <- file.path(
    proc_dir,
    "rasterize_terrain.vrt"
  )
  if (!file.exists(vrt_check) || .force) {
    cli::cli_alert_info("Creating DTM")
    dtm <- lidR::rasterize_terrain(laz_ctg,
                                   res = 1,
                                   algorithm = .alg
    )
  } else {
    cli::cli_alert_info("DTM already exists")
    dtm <- terra::rast(vrt_check)
  }
  return(
    switch(.return_class,
           "SpatRaster" = dtm,
           "character" = terra::sources(dtm)
    )
  )
}
#' @title Generate a normalized point cloud from a laz catalog
#' @description Generate a normalized point cloud from a laz catalog with
#' opinoated settings
#' @param laz_ctg A LAScatalog object
#' @param dtm A SpatRaster object
#' @param .force Force the normalized point cloud to be generated even if it
#' already exists
#' @return A LAScatalog object
lidar_build_norm_pnts <- function(laz_ctg, dtm, .force = FALSE) {
  lid_proc_dir <- dirname(laz_ctg@data$filename[1])
  mod_name <- basename(lid_proc_dir)
  proc_dir <- file.path(lid_proc_dir, paste0(mod_name, "_norm_ctg"))
  lidR::opt_output_files(laz_ctg) <- paste(
    proc_dir,
    paste0(mod_name, "_normal_points_{ID}"),
    sep = "/"
  )
  if (length(list.files(proc_dir)) > 0 || .force) {
    cli::cli_warn(c(
      "!" = "Normalized points already exist. Skipping normalization.",
      "i" = "You should check if the files are present and correct."
    ))
    norm_ctg <- lidR::readLAScatalog(proc_dir)
  } else {
    cli::cli_alert_info("Normalizing points")
    norm_ctg <- lidR::normalize_height(laz_ctg, dtm)
  }
  return(norm_ctg)
}

##Test##
#dtm<- lidar_build_dtm(catalog_list[[2]])
#norm<- lidar_build_norm_pnts(catalog_list[[2]], dtm = dtm)

#### functions for  generating metrics #####


grid_rumple_index <- function(las, res) { # user-defined function
  las <- filter_surfacepoints(las, 1)
  return(pixel_metrics(las, rumple_index(X, Y, Z), res))
}
##Test##
#ctg <- catalog_list[[3]]
#opt_chunk_buffer(ctg) <- 10
#opt_chunk_alignment(ctg) <- c(100, 200)
#options <- list(raster_alignment = 20)
#metrics <- catalog_map(ctg, grid_rumple_index, res = 20, .options = options)
#plot(metrics, col = height.colors(50))



LADCV <- function(z) {
  lad <- try(LAD(z, dz = 1, k = 0.5, z0 = 2))
  if (inherits(lad, "try-error")) {
    return(NA_real_)
  }
  
  if (class(lad) == "numeric") {
    return(NA_real_)
  }
  
  r <- mean(lad$lad, na.rm = TRUE)
  
  return(r)
}


#trial_lad <- pixel_metrics(ctg, ~LADCV(Z), res = 10)



GFCV <- function(z) {
  GF <- try(gap_fraction_profile(z))
  if (inherits(GF, "try-error")) {
    return(NA_real_)
  }
  
  if (class(GF) == "numeric") {
    return(NA_real_)
  }
  
  r <- mean(GF$gf, na.rm = TRUE)
  
  return(r)
}

#trial_GF <- pixel_metrics(ctg, ~GFCV(Z), res = 10)
#terra::plot(trial_GF)

Entropy <- function(las, res) { # user-defined function
  return(pixel_metrics(las, entropy(Z), res))
}

##Test##
#ctg <- catalog_list[[2]]
#opt_chunk_buffer(ctg) <- 10
#opt_chunk_alignment(ctg) <- c(100, 200)
#options <- list(raster_alignment = 20)
#metrics <- catalog_map(ctg, Entropy, res = 20, .options = options)
#plot(metrics, col = height.colors(50))

VCIfunc <- function(las, res) { # user-defined function
  return(pixel_metrics(las, VCI(Z, zmax= max(Z)), res))
}

#ctg <- catalog_list[[2]]
#opt_chunk_buffer(ctg) <- 10
#opt_chunk_alignment(ctg) <- c(100, 200)
#options <- list(raster_alignment = 20)
#metrics <- catalog_map(ctg, VCIfunc, res = 20, .options = options)
#plot(metrics, col = height.colors(50))

# calculate standard metrics
stmetrics <- function(las, res) {
  return(pixel_metrics(las, .stdmetrics, res))
} 
#metrics <- catalog_map(ctg, stmetrics, res = 20, .options = options)
#plot(metrics, col = height.colors(50))



###### function to run all files through all functions #####
run_everything<- function(laz_ctg){
  name<-basename(laz_ctg$filename)
  dtm<- lidar_build_dtm(laz_ctg)
  norm<- lidar_build_norm_pnts(laz_ctg, dtm = dtm)
  opt_chunk_buffer(laz_ctg) <- 5
  opt_chunk_alignment(laz_ctg) <- c(100, 200)
  options <- list(raster_alignment = 20)
  
  lid_proc_dir <- dirname(laz_ctg@data$filename[1])
  mod_name <- basename(lid_proc_dir)
  
  rump_dir <- file.path(lid_proc_dir, paste0(mod_name, "_rumple"))
  lidR::opt_output_files(laz_ctg) <- paste(
    rump_dir,
    paste0(mod_name, "_rumple_{ID}"),
    sep = "/"
  )
  rump_vrt_check <- file.path(
    rump_dir,
    "FUN.vrt"
  )
  if (!file.exists(rump_vrt_check)) {
    cli::cli_alert_info("Creating rumple")
    rump <- catalog_map(laz_ctg, grid_rumple_index, res = 20, .options = options)
  } else {
    cli::cli_alert_info("rumple already exists")
    rump <- terra::rast(rump_vrt_check)
  }
  
  LAD_dir <- file.path(lid_proc_dir, paste0(mod_name, "_LAD"))
  lidR::opt_output_files(laz_ctg) <- paste(
    LAD_dir,
    paste0(mod_name, "_LAD_{ID}"),
    sep = "/"
  )
  lad_vrt_check <- file.path(
    LAD_dir,
    "pixel_metrics.vrt"
  )
  if (!file.exists(lad_vrt_check)) {
    cli::cli_alert_info("Creating LAD")
    lad <- pixel_metrics(laz_ctg, ~LADCV(Z), res = 10)
  } else {
    cli::cli_alert_info("LAD already exists")
    lad <- terra::rast(lad_vrt_check)
  }
  
  
  GF_dir <- file.path(lid_proc_dir, paste0(mod_name, "_GF"))
  lidR::opt_output_files(laz_ctg) <- paste(
    GF_dir,
    paste0(mod_name, "_GF_{ID}"),
    sep = "/"
  )
  GF_vrt_check <- file.path(
    GF_dir,
    "pixel_metrics.vrt"
  )
  if (!file.exists(GF_vrt_check)) {
    cli::cli_alert_info("Creating GF")
    GF <- pixel_metrics(laz_ctg, ~GFCV(Z), res = 10)
  } else {
    cli::cli_alert_info("GF already exists")
    GF <- terra::rast(GF_vrt_check)
  }
  
  ent_dir <- file.path(lid_proc_dir, paste0(mod_name, "_ent"))
  lidR::opt_output_files(laz_ctg) <- paste(
    ent_dir,
    paste0(mod_name, "_ent_{ID}"),
    sep = "/"
  )
  ent_vrt_check <- file.path(
    ent_dir,
    "FUN.vrt"
  )
  if (!file.exists(ent_vrt_check)) {
    cli::cli_alert_info("Creating ent")
    ent <- catalog_map(laz_ctg, Entropy, res = 20, .options = options)
  } else {
    cli::cli_alert_info("ent already exists")
    ent <- terra::rast(ent_vrt_check)
  }
  
  VCI_dir <- file.path(lid_proc_dir, paste0(mod_name, "_VCI"))
  lidR::opt_output_files(laz_ctg) <- paste(
    VCI_dir,
    paste0(mod_name, "_VCI_{ID}"),
    sep = "/"
  )
  VCI_vrt_check <- file.path(
    VCI_dir,
    "FUN.vrt"
  )
  if (!file.exists(VCI_vrt_check)) {
    cli::cli_alert_info("Creating VCI")
    VCI <- catalog_map(laz_ctg, VCIfunc, res = 20, .options = options)
  } else {
    cli::cli_alert_info("VCI already exists")
    VCI <- terra::rast(VCI_vrt_check)
  }
  
  std_dir <- file.path(lid_proc_dir, paste0(mod_name, "_std"))
  lidR::opt_output_files(laz_ctg) <- paste(
    std_dir,
    paste0(mod_name, "_std_{ID}"),
    sep = "/"
  )
  std_vrt_check <- file.path(
    std_dir,
    "FUN.vrt"
  )
  if (!file.exists(std_vrt_check)) {
    cli::cli_alert_info("Creating std")
    std <- catalog_map(laz_ctg, stmetrics, res = 20, .options = options)
  } else {
    cli::cli_alert_info("STD already exists")
    std <- terra::rast(std_vrt_check)
  }
  
  
  comp_dir <- file.path(lid_proc_dir, paste0(mod_name, "_comp.vrt"))
  comp <- rsi::stack_rasters(list(file.path(rump_dir,"FUN.vrt"),
                                  file.path(LAD_dir,"pixel_metrics.vrt"),
                                  file.path(GF_dir,"pixel_metrics.vrt"),
                                  file.path(ent_dir,"FUN.vrt"),
                                  file.path(VCI_dir,"FUN.vrt"),
                                  file.path(std_dir,"FUN.vrt")),
                             band_names = c("rump","LAD", "GF", "ent", "VCI", "zmax", 
                                            "zmean","zsd","zskew","zkurt","zentropy",
                                            "pzabovezmean", "pzabove2","zq5","zq10","zq15", "zq20",        
                                             "zq25" ,"zq30","zq35","zq40","zq45","zq50",       
                                            "zq55","zq60","zq65","zq70","zq75","zq80",        
                                            "zq85","zq90","zq95","zpcum1","zpcum2","zpcum3",      
                                            "zpcum4","zpcum5","zpcum6","zpcum7","zpcum8","zpcum9",      
                                            "itot","imax","imean","isd","iskew","ikurt",      
                                            "ipground","ipcumzq10","ipcumzq30","ipcumzq50","ipcumzq70","ipcumzq90",   
                                            "p1th","p2th","p3th","p4th","p5th","pground",    
                                            "n","area"),
                                 output_filename = paste(comp_dir))
 
  
  return(comp)
} 
  

##Test##
#trial_stack<- run_everything(catalog_list[[2]])
#stack<- terra:: rast(trial_stack)

stack_list<-ctg_list%>%
  map(run_everything)



#### debgging ####