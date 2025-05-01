### A script to use GLMs to model and analyse the mechanistic effects of canopy structure on Soil C

#### set up environment ####
install.packages("readr")
library(performance)
library(sf)
library(dplyr)
library(stringr)
library(stats)
library(tidyr)
library(readr)

##PCA to reduce variables 
# Drop geometry
retained_vars<- read_csv("~/workspace/PhD_work/soil_chapter/Data/variable_presence_count.csv")

## filter for retained vars from ML
fixed_vars <- retained_vars$value
vars_present <- intersect(fixed_vars, colnames(converted_data))

# Select only those columns from the retained in ML
df <- converted_data%>%
  st_drop_geometry()%>%
  select(all_of(vars_present))%>%
  select(-"wmean_percN_5")

# Encode factor columns to numeric
df_numeric <- df %>%
  select(where(is.numeric))%>%
  select(where(~ sd(., na.rm = TRUE) != 0))



# Filter for columns with canopy structure metrics
df_CS <- df_numeric %>%
  select(matches("(_2$|_3$|_4$)"))

#Run PCA on filtered data
rda.out <- vegan::rda(df_numeric, scale = TRUE)
# add scores()
rda_scores <- vegan::scores(rda.out)
# add biplot()
biplot(rda.out, display = c("species"), type = "text")

#Run PCA on CS data
CSrda.out <- vegan::rda(df_CS, scale = TRUE)
# add scores()
CSrda_scores <- vegan::scores(CSrda.out)
# add biplot()
biplot(CSrda.out, display = c("species"), type = "text")

##bi_plot with cos2 scores

bi_plot_func<- function(PCA_scores) {
  
  # Extract species (variable) scores for the first two axes
  species_scores <- as.data.frame(PCA_scores$species[, 1:2])
  colnames(species_scores) <- c("PC1", "PC2")
  species_scores$var <- rownames(species_scores)
  
  # Compute cos2 = (loading^2) / (sum of squares of all loadings for that variable)
  species_scores <- species_scores %>%
    mutate(
      cos2 = PC1^2 + PC2^2
    )
  
  # Optional: scale arrows by cos2 for better visual emphasis
  species_scores <- species_scores %>%
    mutate(
      PC1_scaled = PC1 * cos2,
      PC2_scaled = PC2 * cos2
    )
  
  # Plot biplot with cos2 coloring
  p<-ggplot(species_scores, aes(x = PC1_scaled, y = PC2_scaled, label = var)) +
    geom_segment(aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled, color = cos2),
                 arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
    geom_text(nudge_x = 0.02, nudge_y = 0.02, size = 3) +
    scale_color_viridis_c(option = "C") +
    labs(title = "PCA Biplot with cos² from vegan::rda",
         x = "PC1",
         y = "PC2",
         color = "cos²") +
    theme_minimal()
  return (p)
}

full_plot<- bi_plot_func(rda_scores)
CS_plot<- bi_plot_func(CSrda_scores)

### chatGPT script for using PCA to remove vars with high collinearity
loadings <- as.data.frame(vegan::scores(rda.out, choices = 1:2, display = "species"))
CSloadings <- as.data.frame(vegan::scores(CSrda.out, choices = 1:2, display = "species"))
#Compute cosine similarity matrix (i.e., variable similarity in loading space)
cosine_similarity <- function(x) {
  x <- as.matrix(x)
  sim <- x %*% t(x)
  norms <- sqrt(rowSums(x^2))
  sim / outer(norms, norms)
}
cos_mat <- cosine_similarity(loadings)
CScos_mat <- cosine_similarity(CSloadings)
# 4. Find variable pairs with high similarity (e.g., > 0.9)
highly_similar <- which(abs(cos_mat) > 0.9 & abs(cos_mat) <= 1, arr.ind = TRUE)
CShighly_similar <- which(abs(CScos_mat) > 0.9 & abs(CScos_mat) <= 1, arr.ind = TRUE)
var_pairs <- unique(t(apply(highly_similar, 1, sort)))
CSvar_pairs <- unique(t(apply(CShighly_similar, 1, sort)))
vars_to_remove <- unique(c(rownames(loadings)[var_pairs[, 2]],rownames(CSloadings)[CSvar_pairs[, 2]] ))

# 5. Filter the dataframe to drop redundant variables
df_reduced <- converted_data %>% select(-all_of(vars_to_remove)) %>% st_drop_geometry()

# Output
cat("Removed", length(vars_to_remove), "variables due to high collinearity.\n")
print(vars_to_remove)
