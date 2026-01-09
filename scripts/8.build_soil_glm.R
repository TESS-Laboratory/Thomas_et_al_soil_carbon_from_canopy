## A script to use GLMs to model and analyse the mechanistic effects of canopy structure on Soil C

#### set up environment ####
install.packages("readr")
install.packages("ggrepel")
install.packages("broom")
install.packages("ggeffects")
install.packages("sjPlot")
install.packages("vegan")
install.packages("modelsummary")
install.packages("pandoc")
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
library(modelsummary)

#### Create Plotting theme ####
theme_beautiful <- function() {
  theme_bw(base_size = 13) +
    theme(
      text = element_text(family = "Helvetica"),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 8, color = "black"),
      axis.line.x = element_line(size = 0.3, color = "black"),
      axis.line.y = element_line(size = 0.3, color = "black"),
      axis.ticks = element_line(size = 0.3, color = "black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
      plot.title = element_text(
        size = 8,
        vjust = 1,
        hjust = 0.5,
        color = "black"
      ),
      legend.text = element_text(size = 8, color = "black"),
      legend.title = element_text(size = 8, color = "black"),
      legend.position = c(0.9, 0.9),
      legend.key.size = unit(0.9, "line"),
      legend.background = element_rect(
        color = "black",
        fill = "transparent",
        size = 2,
        linetype = "blank"
      )
    )
}
windowsFonts("Helvetica" = windowsFont("Helvetica")) # Ensure font is mapped correctly

### table set up
set_sf<-function(x){ ifelse(
  x > 1,
  fmt_sprintf("%.2f")(x),
  ifelse(abs(x) <= 1e-3,
         formatC(x, format = "e", digits = 2),
         ifelse(abs(x) <= 1, 
                formatC(x, format = "fg", digits = 2),
                formatC(x, format = "f", digits = 2))
  ))
}
set.seed(42)

##set WD
#online
#fp<- "~/workspace/PhD_work/soil_chapter/Data/"
#local
fp<-"C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data"
#### load data
samples_metrics<- read_rds(file.path(fp, "soil_samples_w_complete_metrics.rds"))

train_data<- samples_metrics|>
  dplyr::rename(percC = percC_5)



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

converted_data <- convert_to_numeric_or_factor(samples_metrics)



#### format table for glm ##### 
df <- converted_data%>%
  st_drop_geometry()%>%
  select(-"percN_5", -"x15N_5", -"fhd_normal_2", -"pground_3", -"...15")|>
  select( -c('VCI_3', 'ent_3', 'zentropy_3'))# , 'GF_3', 'VCI_3'


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

df_numeric <- df_means |>
  mutate(across(where(~!is.numeric(.)), ~ as.numeric(as.factor(.x))))|>  # Convert non-numeric to numeric factors
  select(where(is.numeric))|>                                           # Select all numeric columns
  select(where(~ sd(.x, na.rm = TRUE) != 0))|>
  select(-matches("^zq|cover_z|pai_z|pavd_z|zpcum"))|>
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

### PCA to reduce predictors


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
    if (var == "percC_5") next  # Skip this specific variable
    
    if (!(var %in% selected_vars)) {
      selected_vars <- c(selected_vars, var)
      break
    }
    
    # Stop once we have 20 unique variables
    if (length(selected_vars) >= 20) break
  }}
#df_selected <- df_numeric[, selected_vars]
df_selected <- as.data.frame(st_drop_geometry(df_means))[, selected_vars]

fixed_vars<- colnames(df_selected)
#create a formula for chosen variables
formula_str <- paste("percC_5 ~", paste(fixed_vars, collapse = " + "))
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
long_names<-c("LAD_3" = "Leaf Area Density (Airborne Lidar)",
              "wmean_soc_1"  = "Soil Organic Carbon (Online)",
              "year_of_last_fire_1" = "Year of Most Recent Fire (Online)",
              "imean_3" = "Mean Intensity of Lidar Returns (Airborne Lidar)",
              "p1th_3"  = "Percentage 1st Returns (Airborne Lidar)",
              "wmean_sand_1"  =  "Sand content (Online)",
              "isd_3" = "Standard Deviation of Lidar Return Intensity (Airborne Lidar)",
              "secondary_vegetation_age_2023_1" = "Age of Secondary Vegetation (Online)",
              "classification_2023_1"   = "Classification of Forest Degredation (Online)",
              "itot_3"  = "Total Lidar Return Intensity (Airborne Lidar)",
              "pgap_theta_2" = "Pgap(theta) of Canopy (Spaceborne Lidar)",
              "pzabovezmean_3"  = "Percentage of Returns Above Mean Canopy Height (Airborne Lidar)",
              "imax_3" = "Maximum Intensity of Lidar Returns (Airborne Lidar)",
              "zkurt_3" = "Kurtosis of Canopy Height (Airborne Lidar)",
              "ipground_3" = "Percentage of Intensity Returned by Ground (Airborne Lidar)",
              "clumping.index_4" = "Clumping Index (Hemispherical Photography)",
              "wmean_cec_1" = "Cation Exchange Capacity (Online)",
              "Seasonal_NDVI_Change_1"  = "Seasonall Difference in NDVI (Online)", 
              "canopy.openness_4" = "Canopy Openness (Hemispherical Photography)",
              "wmean_nitrogen_1" = "Total Nitrogen (Online)")
table_s3<-modelsummary(list("20 variable GLM Summary" = glm_model), 
             fmt = set_sf, 
             shape = term ~ model + statistic,
             coef_omit = "Intercept",
             estimate  = "{estimate}",
             statistic = c("{std.error}", "{p.value}"),
             gof_omit = "Num.Obs|BIC|Log.Lik.|F",
             exponentiate = TRUE,
             output = "Plots/Table_S3.docx"
)
TableS3<- tidy(glm_model, exponentiate = TRUE)
write.csv(TableS3, "Plots/TableS3.csv")

# Perform automated backward selection
best_model <- step(glm_model, direction = "backward")

#View the summary of the best model
summary(best_model)
check_model(best_model)
table_s4<-modelsummary(list("Most Parsimonious GLM Summary" = best_model), 
                       fmt = set_sf, 
                       shape = term ~ model + statistic,
                       coef_omit = "Intercept",
                       estimate  = "{estimate}",
                       statistic = c("{std.error}", "{p.value}"),
                       gof_omit = "Num.Obs|BIC|Log.Lik.|F",
                       exponentiate = TRUE,
                       output = "Plots/Table_S4.docx"
)
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

#####  further investigation (should be done manually)

##remove predictors with high colinearity and modelled data 
simple_glm <-glm(formula = percC_5 ~ LAD_3 +  imean_3 + p1th_3 + 
                    itot_3 + pzabovezmean_3 + zkurt_3, 
                 family = Gamma(link = "log"), data = df_means)
check_model(simple_glm)

coef_df <- tidy(
  simple_glm,
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


ggsave("Plots/effect_sizes_simple.png", plot = effect_plot, width = 15, height = 10, dpi = 300)
table_s5<-modelsummary(list("Reduced GLM Summary" = simple_glm), 
                       fmt = set_sf, 
                       shape = term ~ model + statistic,
                       coef_omit = "Intercept",
                       estimate  = "{estimate}",
                       statistic = c("{std.error}", "{p.value}"),
                       gof_omit = "Num.Obs|BIC|Log.Lik.|F",
                       exponentiate = TRUE,
                       coef_rename = long_names,
                       output = "Plots/Table_Ss.docx"
)

TableS5<-tidy(simple_glm, exponentiate = TRUE)
write.csv(TableS5, "Plots/TableS5.csv")

p5<-plot(predict_response(simple_glm, terms= "LAD_3")) +
labs( title = "a)" , x= "Mean Leaf Area Density", y= "Mean %C")+
  theme_minimal(base_size = 24)
p6<-plot(predict_response(simple_glm, terms= "pzabovezmean_3")) +
  labs( title = "b)" , x= "Percentage of Lidar Returns Above Mean Canopy Height", y= "Mean %C")+
  theme_minimal(base_size = 24)
p7<-plot(predict_response(simple_glm, terms= "imean_3")) +
  labs( title = "c)" , x= "Mean Intensity of Lidar Returns", y= "Mean %C")+
  theme_minimal(base_size = 24)
p8<-plot(predict_response(simple_glm, terms= "zkurt_3")) +
  labs( title = "d)" , x= "Kurtosis of Canopy Height Distribution", y= "Mean %C")+
  theme_minimal(base_size = 24)
p9<-plot(predict_response(simple_glm, terms= "p1th_3")) +
  labs( title = "e)" , x= "Percentage 1st Returns", y= "Mean %C")+
  theme_minimal(base_size = 24)
p10<-plot(predict_response(simple_glm, terms= "itot_3 ")) +
  labs( title = "f)" , x= "Total Intensity of Lidar Returns", y= "Mean %C")+
  theme_minimal(base_size = 24)
p<- (p5+p6 )/(p7 +p8)/( p9+ p10)

ggsave('plots/marginal_effects.png', p, width = 20, height = 18, dpi = 300)

