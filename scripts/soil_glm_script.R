#### script to pick out points for analysis and model them 

#### set up environment ####
#install.packages("performance")
#install.packages("GGally")
install.packages("caret")
library(terra)
library(ggplot2)
library(GLMMRR)
library(dplyr)
library(performance)
library(GGally)
library(sf)
library(tidyverse)
library(fastDummies)

#### read in data ####
metrics<- rast("C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Large_Data/combined_metrics_raster.tif")
samples<- read_csv("Data/soil_meta_table.csv")
biomass<-rast("Data/lidar_agb_pred_test_with_zerosV5-0.tif")
rivers<- rast("Data/distance_to_water.tif")
NDVI_var<- rast("C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Large_Data/NDVI_var.tif")
NDVI_grad<- rast("C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Large_Data/NDVI_grad.tif")
#### extract data points ####
# combine observed data
set.ext(biomass, ext(metrics))
set.ext(rivers, ext(metrics))
set.ext(NDVI_grad, ext(metrics))
set.ext(NDVI_var, ext(metrics))

biomass <- resample(biomass, metrics, method = "bilinear")
rivers <- resample(rivers, metrics, method = "bilinear")
NDVI_var <- resample(NDVI_var, metrics, method = "bilinear")
NDVI_grad <- resample(NDVI_grad, metrics, method = "bilinear")

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
  left_join(samples, by = "Codigo") # Use other joins like inner_join if needed

# Print a summary of the merged spatial data frame
print(head(merged_spatial_df))

# Optional: Save the merged spatial data frame to a new shapefile or GeoJSON
#output_path <- "path/to/output/merged_vector.geojson"
#st_write(merged_spatial_df, output_path, delete_dsn = TRUE)


#### full data GLM ####
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
target_column <- "percC"
predictor_columns <- setdiff(names(samples_metrics), target_column)  # Exclude "target"

## remove non numeric variables for pca
numeric_predictors <- sapply(samples_metrics[predictor_columns], is.numeric())  # Check for numeric columns
numeric_columns <- predictor_columns[numeric_predictors]  # Keep only numeric predictors
data_encoded <- dummy_cols(data, select_columns = "", remove_first_dummy = TRUE, remove_selected_columns = TRUE)

# Scale and center the predictor data
scaled_data <- scale(samples_metrics[numeric_columns])

# Step 2: Perform PCA
pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)

# Summary of PCA
summary(pca_result)

# Step 3: Examine Component Contributions
# Add PCA scores and the target column "x" to a new data frame
pca_scores <- as.data.frame(pca_result$x)
pca_scores[[target_column]] <- data[[target_column]]

# Step 4: Analyze Relationship Between Principal Components and x
# Correlation between each PC and the target column "x"
correlations <- sapply(pca_scores[, -ncol(pca_scores)], function(pc) cor(pc, pca_scores[[target_column]], use = "complete.obs"))
print("Correlations between PCs and target column:")
print(correlations)

# Step 5: Visualize PCA Results
# Biplot of PCA
biplot(pca_result, scale = 0)

# Scatter plot of the most influential PCs with "x"
most_influential_pc <- names(which.max(abs(correlations)))  # Identify the most correlated PC
ggplot(pca_scores, aes_string(x = most_influential_pc, y = target_column)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = paste("Scatter Plot of", most_influential_pc, "and", target_column),
       x = most_influential_pc,
       y = target_column) +
  theme_minimal()

# Step 6: Assess Variance Explained
# Proportion of variance explained by each principal component
variance_explained <- summary(pca_result)$importance[2, ]
print("Variance Explained by Each PC:")
print(variance_explained)
#### all RS data GLM ####

#### only open data GLM ####

#### Debugging ####
ggpairs()
