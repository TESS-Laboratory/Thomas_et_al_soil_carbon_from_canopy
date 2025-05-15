## A script to use GLMs to model and analyse the mechanistic effects of canopy structure on Soil C

#### set up environment ####
install.packages("readr")
install.packages("ggrepel")
library(performance)
library(sf)
library(dplyr)
library(stringr)
library(stats)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)

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
  species_scores <- species_scores %>%
    mutate(loading_magnitude = sqrt(PC1_scaled^2 + PC2_scaled^2))
  
  top_vars <- species_scores %>%
    slice_max(loading_magnitude, n = 20)
  
  
  # Plot biplot with cos2 coloring
  p<-ggplot(species_scores, aes(x = PC1_scaled, y = PC2_scaled, label = var)) +
    geom_segment(aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled, color = cos2),
                 arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
    geom_text_repel(nudge_x = 0.02, nudge_y = 0.02, size = 3,max.overlaps = 5, 
                    box.padding = 0.35, point.padding = 0.5) +
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


# Get loadings (rotation matrix: variables x components)
loadings <- as.data.frame(vegan::scores(rda.out, choices = 1:20, display = "species")) # Take absolute values for proximity

# Create an empty vector to collect selected variable names
selected_vars <- c()

# Loop through PCs 1 to 20
for (i in 1:ncol(loadings)) {
  
  # Get sorted variable names by loading (descending) for this PC
  pc<- names(loadings)[1]
  
  ranked_vars <- loadings%>%
    as.data.frame() %>%
    tibble::rownames_to_column("var") %>%
    arrange(desc(.data[[pc]])) %>%
    pull(var)
  # Find the first variable not already selected
  for (var in ranked_vars) {
    if (!(var %in% selected_vars)) {
      selected_vars <- c(selected_vars, var)
      break
    }
  }
  
  # Stop once we have 20 unique variables
  if (length(selected_vars) >= 20) break
}
df_selected <- df[, selected_vars]
df_selected <- as.data.frame(st_drop_geometry(converted_data))[, selected_vars]

fixed_vars<- colnames(df_selected)
#create a formula for chosen variables
formula_str <- paste("wmean_percC_5 ~", paste(fixed_vars, collapse = " + "))
glm_formula <- as.formula(formula_str)

##compare performance of different family and links
glm.gaus.log <- glm(glm_formula, data = converted_data, family= gaussian('log'))
glm.gaus.id <- glm(glm_formula, data = converted_data, family = gaussian('identity'))
glm.gam.log <- glm(glm_formula, data = converted_data, family = Gamma(link="log"))
glm.gam.id<- glm(glm_formula, data = converted_data, family = Gamma(link="identity"))

## review
compare_performance(glm.gam.log, glm.gaus.id, glm.gaus.log, glm.gam.id, verbose = FALSE, rank = TRUE)
summary(glm.gam.log)

check_model(glm.gam.log)

plot(resid(glm.gam.log, type='response'))
lines(resid(glm.gam.log, type='response'), col='red')

##selected formula
glm_model <- glm(glm_formula, data = converted_data, family = Gamma(link="log"))

summary(glm_model)


# Perform automated backward selection
best_model <- step(glm_model, direction = "backward")

#View the summary of the best model
summary(best_model)

## plot effect sizes 
coef_summary <- summary(best_model)$coefficients
conf_int <- confint(best_model)

# Create a data frame for plotting
effects_df <- data.frame(
  Predictor = rownames(coef_summary),
  Estimate = coef_summary[, "Estimate"],
  CI_low = conf_int[, 1],
  CI_high = conf_int[, 2]
)

# remove the intercept to focus only on predictors
effects_df <- effects_df[effects_df$Predictor != "(Intercept)", ]

# Plot
ggplot(effects_df, aes(x = reorder(Predictor, Estimate), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Effect Sizes (Log Scale) from Gamma GLM",
    x = "Predictor",
    y = "Log Coefficient Estimate ± 95% CI"
  )

