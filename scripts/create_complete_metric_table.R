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
library(tidyverse)
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
metrics<- rast(file.path(fp, "combined_metrics_raster.tif"))
#update(x, names=TRUE)
samples<- read_csv(file.path(fp, "soil_meta_table.csv"))
GEDI<- read_csv(file.path(fp, "Gedi_2b_dataframe.csv"))
colnames(GEDI) <- paste0(colnames(GEDI), "_2")
GEDI<-rename(GEDI, Codigo = Codigo_2)

l_metrics<- metrics

# Convert coordinate table to a SpatVector for extraction
coordinates <- vect(samples[, c("plot.x", "plot.y", "Codigo")], geom = c("plot.x", "plot.y"), crs = "EPSG:4326")

#align coordinates to rasters
proj<- crs(metrics)
coordinates <- project(coordinates, proj, partial = TRUE)

# Extract values across all layers in the raster stack
extracted_values <- terra::extract(l_metrics, coordinates, bind = TRUE)

# Convert the SpatVector to an sf object 
spatial_df <- st_as_sf(extracted_values)

# Bind the columns from the original samples table to the spatial data frame
merged_spatial_df <- spatial_df %>%
  left_join(samples, by = "Codigo")%>%
  left_join(GEDI, by = "Codigo")

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
samples_metrics <- clean_headers(merged_spatial_df)%>%
  select(-...1, -Local, -notes, -"Identifier1", -massa, -Ponto, -...15)
samples_metrics <- samples_metrics %>%
  rename_with(~ ifelse(grepl("^\\d", .), paste0("x", .), .))

st_write(samples_metrics, file.path(fp, "soil_samples_w_complete_metrics.csv"), delete_dsn= TRUE)
