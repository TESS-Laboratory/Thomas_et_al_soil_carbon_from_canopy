#### script to pick out points for analysis and model them 

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


#### read in data ####
metrics<- rast("C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Large_Data/combined_metrics_raster.tif")
samples<- read_csv("Data/soil_meta_table.csv")
GEDI<- read_csv("Data/Gedi_2b_dataframe.csv")
colnames(GEDI) <- paste0(colnames(GEDI), "_FA")
GEDI<-rename(GEDI, Codigo = Codigo_FA)
biomass<-rast("Data/lidar_agb_pred_test_with_zerosV5-0.tif")
names(biomass) <- paste0(names(biomass), "_RS")
rivers<- rast("Data/distance_to_water.tif")
names(rivers) <- paste0(names(rivers), "_FA")
NDVI_var<- rast("C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Large_Data/NDVI_var.tif")
names(NDVI_var) <- paste0(names(NDVI_var), "_FA")
NDVI_grad<- rast("C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Large_Data/NDVI_grad.tif")
names(NDVI_grad) <- paste0(names(NDVI_grad), "_FA")
#### extract data points ####
# combine observed data
set.ext(biomass, ext(metrics))
set.ext(rivers, ext(metrics))
set.ext(NDVI_grad, ext(metrics))
set.ext(NDVI_var, ext(metrics))

biomass <- terra::resample(biomass, metrics, method = "bilinear")
rivers <- terra::resample(rivers, metrics, method = "bilinear")
NDVI_var <- terra::resample(NDVI_var, metrics, method = "bilinear")
NDVI_grad <- terra::resample(NDVI_grad, metrics, method = "bilinear")

l_metrics<- c(metrics, biomass, rivers, NDVI_grad, NDVI_var)

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
  #rename_with(~ paste0(., "_FA"), .cols = setdiff(names(spatial_df), "Codigo"))
#st_geometry(merged_spatial_df) <- "geometry_FA_FA"
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
  select(-...1, -Local, -notes, -"Identifier1", -massa, -Ponto, -Age_rectified, -Age.category,-...15)
samples_metrics <- samples_metrics %>%
  rename_with(~ ifelse(grepl("^\\d", .), paste0("x", .), .))

st_write(samples_metrics, "Data/soil_samples_w_complete_metrics.fgb", delete_dsn= TRUE)

####prelim analysis #####
# Create a formula for the GLM where the response column is modeled by all other columns
pairplot <- GGally::ggpairs(samples_metrics, columns = c(2:111, 114, 117:141), cardinality_threshold = 50)
pairplot <- GGally::ggpairs(samples_metrics, columns = c(132:133,136:140, 120), cardinality_threshold = 50)
primary_var <- "percC"
pvar_pos <- match(primary_var, pairplot$xAxisLabels)
plots <- lapply(1:pairplot$ncol, function(j){ getPlot(pairplot, i = pvar_pos, j = j)})
ggmatrix(
  plots,
  nrow = 10,
  ncol = 10,
  xAxisLabels = pairplot$xAxisLabels,
  yAxisLabels = primary_var
)
print(pairplot)
ggsave("var132_140.png", pairplot, path = "Data/Outputs/Pairwaise_plots", width = 35, height = 30, units = "cm")
getPlot(pairplot, i= 3, j= 1)
p_ <- GGally::print_if_interactive

# Rearrange to put column 'x' first
df_rearranged <- samples_metrics[, c("X.C", setdiff(names(samples_metrics), "X.C"))]
gplot <- GGally::ggpairs(df_rearranged, columns= 1)
gplot$nrow <- 1
gplot$yAxisLabels <- df_rearranged$yAxisLabels[1]
print(gplot)
# Create a formula for the GLM where the response column is modeled by all other columns
# Filter out columns with only one unique level (constant columns)
filtered_table <- samples_metrics[, sapply(samples_metrics, function(col) length(unique(col)) > 1)]
formula <- as.formula(paste(response_column, "~ ."))
glm_model <- glm(X.C~ actual.LAI+ `bdod_0-5cm_mean`, data = filtered_table, family = gaussian)
FD_glm<- glm()


#### PCA Analysis ####

## subset columns for numerical variables

metrics_num <- samples_metrics[,c("bdod_0.5cm_mean","cec_0.5cm_mean",
                        "clay_0.5cm_mean","nitrogen_0.5cm_mean","ocs_0.30cm_mean" ,
                        "phh2o_0.5cm_mean",
                        "sand_0.5cm_mean","soc_0.5cm_mean", "rump", "LAD","GF",
                        "ent","VCI",  "zmean", "HydroRiver_raster", "NDVI_variance",
                        "15N","percN","13C", "percC", "CperN", "effective.LAI",          
                        "actual.LAI","clumping.index",         
                        "LXG1", "MTA", "canopy.openness", "min_distance", "Degradation", "Profundidade_cm_")]

metrics_num$Degradation<- as.factor(metrics_num$Degradation)
metrics_num$Profundidade_cm_<- as.factor(metrics_num$Profundidade_cm_)
metrics_num<-data.frame(metrics_num)
metrics_num<- na.omit(metrics_num)

## filter for deapths we are interested in
Fmetrics_num<- filter(metrics_num, Profundidade_cm_ == "0-5" | Profundidade_cm_ =="5-10"| Profundidade_cm_== "10-20" | Profundidade_cm_ == "20-30"  )

## plot without categorical data 
rda.out <- vegan::rda(Fmetrics_num[,-c(29,30,31)], scale = TRUE)
# add scores()
rda_scores <- scores(rda.out)
# add biplot()
biplot(rda.out, type = "text")
## group by depth ### this doesnt work yet
ordihull(rda.out,
         group = Fmetrics_num$Profundidade_cm_,
         col = 1:11,
         lty = 1:11,
         lwd = c(3,6), 
         label = TRUE)


#####  mlr3 ecosystem trial #####

##convert datatable to a task (must remove geometry)
non_sf<-samples_metrics%>%
  st_drop_geometry()%>%
  select(!"Codigo")

## make everything numeric or factor
convert_to_numeric_or_factor <- function(df) {
  df %>%
    mutate(across(everything(), ~ {
      if (suppressWarnings(all(!is.na(as.numeric(.))))) {
        as.numeric(.)
      } else {
        as.factor(.)
      }
    }))
}

# Apply the conversion function
converted_data <- convert_to_numeric_or_factor(non_sf)

T_metrics= as_task_regr(converted_data, target = "percC",
                         id = "soils")
mlr3viz::autoplot(T_metrics, type = "pairs")

## training model
# create a regression tree learner
#any seed will do. just for consistency
set.seed(42)
#ranger is Random forests
trial_learner<-lrn("regr.ranger")
#dummy learner
lrn_featureless = lrn("regr.featureless")

#feature selection
flt_gain = flt("information_gain")
flt_gain$calculate(T_metrics)
variables_by_information_gain<-as.data.table(flt_gain)

T_metrics$select(c("LAD", "clay_0.5cm_mean", "ent",
                 "nitrogen_0.5cm_mean"))
instance = fselect(
  fselector = fs("sequential"),
  task =  T_metrics,
  learner = trial_learner,
  resampling = rsmp("cv", folds = 3),
  measure = msr("regr.mae")
)
dt = as.data.table(instance$archive)
dt[batch_nr == 1, 1:5]## in this example nitrogen is the top variable
autoplot(instance, type = "performance") ## all 4 contribute something but best is 2 levels
dt[batch_nr == 2, 1:5] ##LAD is 2nd most important in this example 
instance$result_feature_set ##alphabetical not best performing, but only prints the best ones (which in this case is all)

## split data to train and test-- holdout method
splits = partition(T_metrics)
#train
lrn_featureless$train(T_metrics, splits$train)
trial_learner$train(T_metrics, splits$train)
trial_learner$train(T_metrics, row_ids = splits$train)
#predict
prediction = trial_learner$predict(T_metrics, row_ids = splits$test)
prediction
autoplot(prediction)

##measure model accuracy with mean absolute error
measures = msrs(c("regr.mse", "regr.mae"))
trial_learner$predict(T_metrics, row_ids = splits$test)$score(measures)
lrn_featureless$predict(T_metrics, row_ids = splits$test)$score(measures)


##https://mlr3book.mlr-org.com/chapters/chapter2/data_and_basic_modeling.html
##start here with 2.6 for column task roles 
##3.2.1 Constructing a Resampling Strategy
##3.3.1 BENCHMARKING!! important
