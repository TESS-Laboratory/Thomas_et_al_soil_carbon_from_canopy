#### script to combine all data sets into one sf table for modelling

#### set up environment ####
#install.packages("performance")
#install.packages("GGally")
#install.packages("fastDummies")
#install.packages("vegan")
#install.packages("mlr3")
#install.packages("mlr3viz")
#install.packages("mlr3learners")
#install.packages("ranger")
#install.packages("mlr3filters")
#install.packages("FSelectorRcpp")
#install.packages("mlr3fselect")
#install.packages("mlr3pipelines")
#install.packages("praznik")
#install.packages("mlr3tuningspaces")
library(scales)
library(terra)
library(ggplot2)
library(GLMMRR)
library(dplyr)
library(performance)
library(GGally)
library(sf)
library(tidyverse)#used
library(fastDummies)
library(vegan)
library(mlr3)
library(mlr3viz)
library(mlr3learners)
library(ranger)
library(mlr3filters)
library(FSelectorRcpp)
library(mlr3fselect)
library(mlr3filters)
library(mlr3pipelines)
library(praznik)
library(mlr3tuningspaces)

###set WD
#workspace

#local
fp<- "C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data"

#### read in data ####
l_metrics<- rast(file.path(fp, "combined_metrics_raster.tif"))

samples<- read_csv(file.path(fp, "soil_meta_table.csv"))
GEDI<- read_csv(file.path(fp, "Gedi_2b_dataframe.csv"))
colnames(GEDI) <- paste0(colnames(GEDI), "_2")
GEDI<-rename(GEDI, Codigo = Codigo_2)

## summarise soil samples 

filter_values <- c("0-5","5-10","10-20","20-30")
filtered_data <- samples |>
  dplyr::mutate(`Profundidade (cm)_5` = as.character(`Profundidade (cm)_5`))|>
  dplyr::filter( `Profundidade (cm)_5` %in% filter_values)
## add weights for means
weight_table <- filtered_data %>%
  dplyr::mutate(weight = dplyr::case_when(
    `Profundidade (cm)_5` == "0-5" ~ 1,
    `Profundidade (cm)_5` =="5-10" ~ 1,
    `Profundidade (cm)_5` == "10-20" ~ 2,
    `Profundidade (cm)_5` == "20-30" ~2
  ))

# Summarize the data by unique IDs with a weighted mean of SOC
summary_samples <- weight_table %>%
  dplyr::group_by(id_clean_5) %>%
  dplyr::summarise(across(everything(), first),
                   across(where(is.numeric), ~ weighted.mean(.x, weight))
                   )|>
  select(-c(`Profundidade (cm)_5`, "id_clean_5", ...1, Local_5, notes_5, `Identifier 1_5`, massa_5, Ponto_5))


# Convert coordinate table to a SpatVector for extraction
coordinates <- vect(summary_samples[, c("plot.x_4", "plot.y_4", "Codigo_5")], geom = c("plot.x_4", "plot.y_4"), crs = "EPSG:4326")

#align coordinates to rasters
proj<- crs(l_metrics)
coordinates <- project(coordinates, proj, partial = TRUE)

# Extract values across all layers in the raster stack
extracted_values <- terra::extract(l_metrics, coordinates, bind = TRUE)

# Convert the SpatVector to an sf object 
spatial_df <- st_as_sf(extracted_values)

# Bind the columns from the original samples table to the spatial data frame
merged_spatial_df <- spatial_df %>%
  left_join(summary_samples, by = "Codigo_5")%>%
  left_join(GEDI, by = c("Codigo_5" = "Codigo"))

# clean table headers
clean_headers <- function(df) {
  clean_names <- gsub("[[:space:]]", "", colnames(df))  # Remove all spaces
  clean_names <- gsub("\\(|\\)", "_", clean_names)  # Remove all brackets
  clean_names<- gsub("%", "perc", clean_names) #remove percent sign
  clean_names<- gsub("/", "per", clean_names)
  clean_names<- gsub("²", "sq", clean_names)
  colnames(df) <- clean_names
  return(df)
}

# Apply the function to the sample data frame
samples_metrics <- clean_headers(merged_spatial_df)
samples_metrics <- samples_metrics %>%
  rename_with(~ ifelse(grepl("^\\d", .), paste0("x", .), .))

train_data_clean<- dplyr::select(samples_metrics, -c("rv_2","shot_number_2","plot.x_4", 
                                                     "plot.y_4",  "scale_1" ,
                                                     ,"date_time_2" ,"light.conditions_4","n_3",
                                                     "GEDI.footprints_4" ,"CperN_5","x13C_5",
                                                     "min_distance_4","Codigo_5", "Age.category_4" ,
                                                     "Age_4", "Age_rectified_4" ,"State_4" ,"Degradation_4"))|>
  dplyr::rename(min_distance_4 = min_distance_1)|>
  filter(st_geometry_type(geometry) == "POINT")



write_rds(train_data_clean, file.path(fp, "soil_samples_w_complete_metrics.rds"))
