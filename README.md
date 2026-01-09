## predicting soil carbon from canopy structure metrics
This repository contains summary data and the analytical codebase of manuscript: 
# Soil carbon predictions across the landscape using remotely- sensed canopy structure measurements in southern Amazonia 

### contributors: Jessica P. Thomas, Andrew M. Cunliffe, Hugh A. Graham, Tom Powell, Plínio B. Camargo, Ted R. Feldpausch.

currently unpublished

Data available at: [(https://doi.org/10.5281/zenodo.17944261)](https://doi.org/10.5281/zenodo.17944261)

To investigate how soil carbon can be predicted using canopy structure metrics, this study addressed the following research questions:

(1) How does the predictive accuracy of  soil carbon models change with different levels of data availability?  
(2) What can we infer about the mechanistic relationships between key predictors and soil carbon?

We selected 142 locations for soil sampling (Figure 1), that were within the footprints sampled by the Global Ecosystem Dynamics Investigation (GEDI) sensor and within the airborne lidar coverage. These locations were stratified to sample different forest ages, degradation histories, distance from rivers, and distance to forest edge.

## This repo contains the following folders: 
Scripts- where all the scripts used in R to collect, clean and analyse data

Data- where the raw and output data used in the analysis were stored

Plots- where any plots resulting from the analysis were saved. 

## The scripts should follow a sequence of:
### 1) LAI_images_to_dataframe.R 
This script processes images taken with a hemispherical camera lens and uses the package hemispheR to convert them into dataframe of coordinates with the relative LAI, canopy openness, an

### 2) lidar_analysis.R
Lidar analysis script to get canopy structure metrics from airborne lidar data using lidR package  (best run on a workstation)

### 3) identify_rivers_raster.R
Script to identify rivers and create the nearest distance to water raster

### 4) compile_rasters.R
script to stack raster layers from GEE, QGIS, lidar_analysis and indentify_rivers_raster into one object

### 5) create_complete_metric_table.R
Script to extract data from the raster stack and combine all data sets into one sf table for modelling with different lables for variables from different sources

### 6) mlr3_modular.R
Script to automate the process of taking each combination of variables and automatically run them through a simultaneous autotuner and feature selection. Then this script runs a comparison of the best models for each variable combination. 

### 7) v2_soil_glm.R
This script is used for the mixed effects modelling to make an inference of the mechanistic relationships between canopy structure and soil Carbon
