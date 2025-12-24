## A script to use GLMs to model and analyse the mechanistic effects of canopy structure on Soil C

#### set up environment ####
install.packages("readr")
install.packages("ggrepel")
install.packages("broom")
install.packages("ggeffects")
install.packages("sjPlot")
install.packages("vegan")
install.packages("gtsummary")
library(gtsummary)
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
library(vegan)
library(ggeffects)
set.seed(42)

##set WD
#online
#fp<- "~/workspace/PhD_work/soil_chapter/Data/"
#local
fp<-"C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data"
#### load data
samples_metrics<- read_sf(file.path(fp, "soil_samples_w_complete_metrics_2.fgb"))
samples_metrics <- samples_metrics|>
  rename_with(~ sub("^wmean_", "", .x))
train_data_clean<- dplyr::select(samples_metrics, -c("acd_lidar_1", "rv_2","shot_number_2","plot.x", 
                                                     "plot.y",  "scale_1" ,
                                              ,"date_time_2" ,"light.conditions","n_3","GEDI.footprints" ,"CperN","x13C","lyr.1_1","Codigo", "Age.category" ,"Age", "Age_rectified" ,"State" ,"Degradation"))  


train_data<- train_data_clean|>
  dplyr::rename(effective.LAI_4 = effective.LAI,
                actual.LAI_4 = actual.LAI,
                clumping.index_4 = clumping.index,
                LXG1_4 = LXG1,
                LXG2_4 = LXG2,
                MTA_4 = MTA,
                canopy.openness_4 = canopy.openness,
                min_distance_4= min_distance,
                Profundidade_cm_5 =  Profundidade_cm_,
                x15N_5 = x15N,
                percN_5 = percN)


filter_values <- c("0-5","5-10","10-20","20-30")
filtered_data <- train_data |>
  dplyr::mutate(Profundidade_cm_5 = as.character(Profundidade_cm_5))|>
  dplyr::filter( Profundidade_cm_5 %in% filter_values)
## add weights for means
weight_table <- filtered_data %>%
  dplyr::mutate(weight = dplyr::case_when(
    Profundidade_cm_5== "0-5" ~ 1,
    Profundidade_cm_5 =="5-10" ~ 1,
    Profundidade_cm_5 == "10-20" ~ 2,
    Profundidade_cm_5 == "20-30" ~2
  ))

# Summarize the data by unique IDs with a weighted mean of SOC
summary_data <- weight_table %>%
  dplyr::group_by(id_clean) %>%
  dplyr::summarise( percC= sum(percC * weight) / sum(weight),
                    across(
                      .cols = -c(percC, weight, geometry),
                      .fns  = first
                    ),
                    .groups = "drop"
  )|>
  select(-c("Profundidade_cm_5", "id_clean"))


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

converted_data <- convert_to_numeric_or_factor(summary_data)

##make sure headers work in mlr3 ecosystem
colnames(converted_data)<-make.names(colnames(converted_data))



#### format table for glm ##### 
df <- converted_data%>%
  st_drop_geometry()%>%
  select(-"percN_5", -"x15N_5", -"fhd_normal_2", -"pground_3")

# weighted means of soil grids
compute_rowwise_weighted_means <- function(df) {
  filter_values<- c("0.5cm", "5.15cm", "15.30cm")
  
  df_long <- df %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(
      cols = matches("*_mean_1$"),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      value_type = str_extract(variable, "[^_]+"),
      depth_range = str_extract(variable, "\\d+\\.\\d+cm"),
      depth_top = as.numeric(str_extract(depth_range, "^\\d+")),
      depth_bottom = as.numeric(str_extract(depth_range, "(?<=\\.)\\d+(?=cm)")),
      depth_thickness = depth_bottom - depth_top
    )|>
      dplyr::filter(depth_range %in% filter_values)
  
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
    select(-matches("*_mean_1$")) %>%
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
  select(where(~ sd(., na.rm = TRUE) != 0))|>
  select(-matches("^zq|cover_z|pai_z|pavd_z"))|>
  drop_na()



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
    geom_text_repel(nudge_x = 0.02, nudge_y = 0.02, size = 8,max.overlaps = 5, 
                    box.padding = 0.35, point.padding = 0.5) +
    scale_color_viridis_c(option = "C") +
    labs(title = "B",
         x = "PC1",
         y = "PC2",
         color = "cos²") +
    theme_minimal()
  return (p)
}

full_plot<- bi_plot_func(rda_scores)
ggsave("Plots/PCA_full_plot.png", plot = full_plot, width = 20, height = 18, dpi = 300)
CS_plot<- bi_plot_func(CSrda_scores)
ggsave("Plots/PCA_CS_plot.png", plot = CS_plot, width = 20, height = 18, dpi = 300)

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
    if (var == "percC") next  # Skip this specific variable
    
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
formula_str <- paste("percC ~", paste(fixed_vars, collapse = " + "))
glm_formula <- as.formula(formula_str)

##compare performance of different family and links
glm.gaus.log <- glm(glm_formula, data = df_means, family= gaussian('log'))
glm.gaus.id <- glm(glm_formula, data = df_means, family = gaussian('identity'))
glm.gam.log <- glm(glm_formula, data = df_means, family = Gamma(link="log"))
glm.gam.id<- glm(glm_formula, data = df_means, family = Gamma(link="identity"))

## review
compare_performance(glm.gam.log, glm.gaus.id, glm.gaus.log, glm.gam.id, verbose = FALSE, rank = TRUE)
summary(glm.gam.log)

check_model(glm.gam.log)

plot(resid(glm.gam.id, type='response'))
lines(resid(glm.gam.id, type='response'), col='red')

##selected formula
glm_model <- glm(glm_formula, data = df_means, family = Gamma(link="log"))

summary(glm_model)
TableS3<- tidy(glm_model, exponentiate = TRUE)
write.csv(TableS3, "Plots/TableS3.csv")

# Perform automated backward selection
best_model <- step(glm_model, direction = "backward")

#View the summary of the best model
summary(best_model)
TableS4<- tidy(best_model, exponentiate = TRUE)
write.csv(TableS4, "Plots/TableS4.csv")

## plot effect sizes 
coef_df <- tidy(
  best_model,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  filter(term != "(Intercept)")

# Plot

effect_plot<-ggplot(coef_df, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2
  ) +
  scale_x_log10() +
  labs(
    x = "Multiplicative effect on percC (log scale)",
    y = NULL,
    title = "Effect sizes from Gamma(log) model",
    subtitle = "Points show exp(coef); bars are 95% confidence intervals"
  ) +
  theme_minimal(base_size = 24)


ggsave("Plots/effect_sizes.png", plot = effect_plot, width = 15, height = 10, dpi = 300)

#####  further investigation

simple_glm <- glm(percC~ min_distance_4+ isd_3 +zkurt_3+imean_3+p1th_3+Seasonal_NDVI_Change_1   ,data = df_means, family = Gamma(link="log"))


TableS5<-tidy(simple_glm, exponentiate = TRUE)
write.csv(TableS5, "Plots/TableS5.csv")

p5<-plot(predict_response(simple_glm, terms= "LAD_3")) +
labs( title = "a)" , x= "Mean Leaf Area Density", y= "Mean %C")
p6<-plot(predict_response(simple_glm, terms= "min_distance_4")) +
  labs( title = "b)" , x= "Minimum Distance from Forest Edge (m)", y= "Mean %C")
p7<-plot(predict_response(simple_glm, terms= "isd_3")) +
  labs( title = "c)" , x= "Standard Deviation of Intensity", y= "Mean %C")
p8<-plot(predict_response(simple_glm1, terms= "zskew_3")) +
  labs( title = "d)" , x= "Skew of Canopy Height", y= "Mean %C")
p<- p5+p6 +p7 +p8

ggsave('plots/model_inference_gradients.png', p, width = 20, height = 18, dpi = 300)

