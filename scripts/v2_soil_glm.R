## A script to use GLMs to model and analyse the mechanistic effects of canopy structure on Soil C

#### set up environment ####
install.packages("readr")
install.packages("ggrepel")
install.packages("broom")
install.packages("ggeffects")
install.packages("sjPlot")
library(performance)
library(sf)
library(dplyr)
library(stringr)
library(stats)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(purrr)
library(broom)
library(ggeffects)
set.seed(42)

##set WD
#online
#fp<- "~/workspace/PhD_work/soil_chapter/Data/"
#local
fp<-"C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data/Data"
#### load data
samples_metrics<- read_sf(file.path(fp, "soil_samples_w_complete_metrics.fgb"))
samples_metrics$wmean_min_distance_4<- samples_metrics$wmean_min_distance_4*1000
#train_data_rm<- select(train_data, where(~!any(is.na(.))))
train_data_clean<- select(samples_metrics, -"wmean_acd_lidar_3", -"wmean_GF_3")%>%
  drop_na()

convert_to_numeric_or_factor <- function(df) {
  df %>%
    mutate(across(where(~ !inherits(., "sfc")), ~ {  # Ignore geometry columns
      num_col <- suppressWarnings(as.numeric(.))
      if (all(is.na(.) | !is.na(num_col))) {
        num_col
      } else {
        as.factor(.)
      }
    }))
}

converted_data <- convert_to_numeric_or_factor(train_data_clean)

##make sure headers work in mlr3 ecosystem
colnames(converted_data)<-make.names(colnames(converted_data))

retained_vars<- read_csv("~/workspace/PhD_work/soil_chapter/Data/variable_presence_count.csv")

#### format table for glm ##### 
df <- converted_data%>%
  st_drop_geometry()%>%
  select(-"wmean_percN_5", -"wmean_x15N_5", -"wmean_fhd_normal_2", -"wmean_scale_1", -"wmean_shot_number_2", -"wmean_rv_2", -"wmean_area_3", -"wmean_pground_3", -"wmean_plot.x", -"wmean_plot.y", -"wmean_x13C_5")

# weighted means of soil grids
compute_rowwise_weighted_means <- function(df) {
  df_long <- df %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(
      cols = matches("^wmean_.*_mean_1$"),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      value_type = str_extract(variable, "(?<=wmean_)[^_]+"),
      depth_range = str_extract(variable, "\\d+\\.\\d+cm"),
      depth_top = as.numeric(str_extract(depth_range, "^\\d+")),
      depth_bottom = as.numeric(str_extract(depth_range, "(?<=\\.)\\d+")),
      depth_thickness = depth_bottom - depth_top
    )
  
  df_weighted <- df_long %>%
    group_by(row_id, value_type) %>%
    summarise(
      weighted_mean = weighted.mean(value, depth_thickness, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = value_type,
      values_from = weighted_mean,
      names_glue = "wmean_{value_type}_1"
    ) %>%
    arrange(row_id)
  
  df_base <- df %>%
    select(-matches("^wmean_.*_mean_1$")) %>%
    mutate(row_id = row_number())
  
  result <- df_base %>%
    left_join(df_weighted, by = "row_id") %>%
    select(-row_id)
  
  return(result)
}

df_means <- compute_rowwise_weighted_means(df)

df_numeric <- df_means %>%
  mutate(across(where(~ !is.numeric(.)), ~ as.numeric(as.factor(.)))) %>%  # Convert non-numeric to numeric factors
  select(where(is.numeric)) %>%                                            # Select all numeric columns
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
    if (var == "wmean_percC_5") next  # Skip this specific variable
    
    if (!(var %in% selected_vars)) {
      selected_vars <- c(selected_vars, var)
      break
    }
    
    # Stop once we have 20 unique variables
    if (length(selected_vars) >= 20) break
  }}
df_selected <- df_numeric[, selected_vars]
df_selected <- as.data.frame(st_drop_geometry(df_means))[, selected_vars]

fixed_vars<- colnames(df_selected)
#create a formula for chosen variables
formula_str <- paste("wmean_percC_5 ~", paste(fixed_vars, collapse = " + "))
glm_formula <- as.formula(formula_str)

##compare performance of different family and links
glm.gaus.log <- glm(glm_formula, data = df_means, family= gaussian('log'))
glm.gaus.id <- glm(glm_formula, data = df_means, family = gaussian('identity'))
glm.gam.log <- glm(glm_formula, data = df_means, family = Gamma(link="log"))
glm.gam.id<- glm(glm_formula, data = df_means, family = Gamma(link="identity"))

## review
compare_performance(glm.gam.log, glm.gaus.id, glm.gaus.log, glm.gam.id, verbose = FALSE, rank = TRUE)
summary(glm.gam.id)

check_model(glm.gam.id)

plot(resid(glm.gam.id, type='response'))
lines(resid(glm.gam.id, type='response'), col='red')

##selected formula
glm_model <- glm(glm_formula, data = df_means, family = Gamma(link="identity"))

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
effects_df <- effects_df %>%
  mutate(
    Estimate = sign(Estimate) * sqrt(abs(Estimate)),
    CI_low = sign(CI_low) *sqrt(abs(CI_low)),
    CI_high = sign(CI_high) *sqrt(abs(CI_high))
  )

# remove the intercept to focus only on predictors
effects_df <- effects_df[effects_df$Predictor != "(Intercept)", ]


# Plot
effect_plot<- ggplot(effects_df, aes(x = reorder(Predictor, Estimate), y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Effect Sizes from GLM",
    x = "Predictor",
    y = "Coefficient Estimate ± 95% CI (sqrt) "
  )+
  theme(text=element_text(size=16))
tidy(best_model, exponentiate = TRUE)


ggsave("Plots/effect_sizes.png", plot = effect_plot, width = 20, height = 18, dpi = 300)

#####  further investigation

simple_glm1 <- glm(wmean_percC_5~ wmean_min_distance_4* wmean_isd_3 *wmean_LAD_3 +wmean_zskew_3 ,data = df_means, family = Gamma(link="identity"))
simple_glm2 <- glm(wmean_percC_5~ wmean_min_distance_4*wmean_isd_3 *wmean_LAD_3 *wmean_zskew_3 ,data = df_means, family = Gamma(link="identity"))
anova(simple_glm1, simple_glm2, test = "Chisq")
summary(simple_glm2)
tidy(simple_glm2, exponentiate = TRUE)
check_model(simple_glm)

p1<-ggplot(df_means)+
  aes(x= wmean_min_distance_4, y= wmean_percC_5)+
  geom_point(alpha= 0.3)+
  geom_smooth(method= "glm", method.args = list(family=gaussian(link="log")))+
  theme_beautiful()
p2<-ggplot(df_means)+
  aes(x= wmean_isd_3, y= wmean_percC_5)+
  geom_point(alpha= 0.3)+
  geom_smooth(method= "glm", method.args = list(family=gaussian(link="log")))+
  theme_beautiful()
p3<-ggplot(df_means)+
  aes(x=  wmean_LAD_3, y= wmean_percC_5)+
  geom_point(alpha= 0.3)+
  geom_smooth(method= "glm", method.args = list(family=gaussian(link="log")))+
  theme_beautiful()
p4<-ggplot(df_means)+
  aes(x=  wmean_zskew_3, y= wmean_percC_5)+
  geom_point(alpha= 0.3)+
  geom_smooth(method= "glm", method.args = list(family=gaussian(link="log")))+
  theme_beautiful()

(p1+p2)/(p3 +p4)


p5<-plot(predict_response(simple_glm, terms= "wmean_LAD_3")) +
labs( title = NULL , x= "Mean Leaf Area Density", y= "Mean %C")
p6<-plot(predict_response(simple_glm, terms= "wmean_min_distance_4")) +
  labs( title = NULL , x= "Minimum Distance from Forest Edge (m)", y= "Mean %C")
p7<-plot(predict_response(simple_glm, terms= "wmean_isd_3")) +
  labs( title = NULL , x= "Standard Deviation of Intensity", y= "Mean %C")
p8<-plot(predict_response(simple_glm, terms= "wmean_zskew_3")) +
  labs( title = NULL , x= "Skew of Canopy Height", y= "Mean %C")
p5+p6 +p7 +p8
