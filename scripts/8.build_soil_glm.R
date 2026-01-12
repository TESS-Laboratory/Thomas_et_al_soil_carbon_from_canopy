## A script to use GLMs to model and analyse the mechanistic effects of canopy structure on Soil C

#### set up environment ####
install.packages("readr")
install.packages("ggrepel")
install.packages("broom")
install.packages("ggeffects")
install.packages("sjPlot")
install.packages("modelsummary")
install.packages("pandoc")
#install.packages("factoextra")
library(factoextra)
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
      legend.text = element_text(size = 18, color = "black"),
      legend.title = element_text(size = 18, color = "black"),
      legend.position = c(0.9, 0.9),
      legend.key.size = unit(3, "line"),
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
  select(-"percN_5", -"x15N_5",  -"...15", -"area_3", - "acd_lidar_3")



# weighted means of soil grids
compute_rowwise_weighted_means <- function(df) {
  
  df_long <- df |>
    mutate(row_id = row_number()) |>
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
    dplyr::filter(depth_bottom<31)
  
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
  select(-matches("^zq|cover_z|pai_z|pavd_z|zpcum|iskew|ikurt|ipcum|itot"))|>
  drop_na()



# Filter for columns with canopy structure metrics
df_CS <- df_numeric %>%
  select(matches("(_2$|_3$|_4$)"))

#Run PCA on filtered data
pca <- prcomp(df_numeric, center = TRUE, scale. = TRUE)

pca_CS<-prcomp(df_CS, center = TRUE, scale. = TRUE)

# Plot biplot with contribution coloring
full_plot<-fviz_pca_var(
  pca,
  repel = TRUE,
  geom = c("arrow", "text"),
  col.var =  "contrib",
  label = "var"
)+
  scale_color_viridis_c(option = "C")+
  labs(title = "A)", color = "Contribution")+
  theme(legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.position = "right",
        legend.key.size = unit(3, "line"),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        plot.title = element_text(
          size = 24,
          color = "black"
        ))
full_plot$layers[[3]]$aes_params$alpha <- 0
full_plot$layers[[2]]$aes_params$linewidth <- 1
full_plot$layers[[1]]$aes_params <- list( nudge_x = 0.02, nudge_y = 0.02, size = 8,max.overlaps = 5, 
                                          box.padding = 0.35, point.padding = 0.5, colour = "black")
full_plot$layers<- c(full_plot$layers[[2]], full_plot$layers[[4]], full_plot$layers[[5]], full_plot$layers[[1]])

ggsave("Plots/PCA_full_plot.png", plot = full_plot, width = 20, height = 20, dpi = 300)

CS_plot<-fviz_pca_var(
  pca_CS,
  repel = TRUE,
  geom = c("arrow", "text"),
  col.var =  "contrib",
  label = "var"
)+
  scale_color_viridis_c(option = "C")+
  labs(title = "B)", color = "Contribution")+
  theme(legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.position = "right",
        legend.key.size = unit(3, "line"),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        plot.title = element_text(
          size = 24,
          color = "black"
        ))
CS_plot$layers[[3]]$aes_params$alpha <- 0
CS_plot$layers[[2]]$aes_params$linewidth <- 1
CS_plot$layers[[1]]$aes_params <- list( nudge_x = 0.02, nudge_y = 0.02, size = 8,max.overlaps = 5, 
                                        box.padding = 0.35, point.padding = 0.5, colour = "black")
CS_plot$layers<- c(CS_plot$layers[[2]], CS_plot$layers[[4]], CS_plot$layers[[5]], CS_plot$layers[[1]])

ggsave("Plots/PCA_CS_plot.png", plot = CS_plot, width = 20, height = 18, dpi = 300)




### PCA to reduce predictors

loadings <- as.data.frame(pca$rotation)
loadings_abs<- abs(loadings)|>
  tibble::rownames_to_column("var")



selected_vars <- c()

# Loop through PCs 1 to 20
for (i in 2:ncol(loadings_abs)) {
  
  # Get sorted variable names by loading (descending) for this PC
  pc<- names(loadings_abs[i])
  
  ranked_vars <- loadings_abs%>%
    arrange(desc(.data[[pc]])) %>%
    pull(var)
  # Find the first variable not already selected
  for (var in ranked_vars) {
    if (var == "percC_5") next  # Skip this specific variable
    
    if (!(var %in% selected_vars)) {
      selected_vars <- c(selected_vars, var)
      break
    }
    
  }# Stop once we have 20 unique variables
  if (length(selected_vars) >= 20) break
}
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
simple_glm <-glm(formula = percC_5 ~ zsd_3 + zskew_3 + zkurt_3 + 
                   isd_3  + p1th_3,
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
                       output = "Plots/Table_Ss.docx"
)

TableS5<-tidy(simple_glm, exponentiate = TRUE)
write.csv(TableS5, "Plots/TableS5.csv")

p5<-plot(predict_response(best_model, terms= "zsd_3")) +
  labs( title = "a)" , x= "Standard Deviation of Canopy Height", y= "Mean %C")+
  theme_minimal(base_size = 24)
p6<-plot(predict_response(best_model, terms= "zskew_3")) +
  labs( title = "b)" , x= "Skew of Canopy Height", y= "Mean %C")+
  theme_minimal(base_size = 24)
p7<-plot(predict_response(best_model, terms= "isd_3")) +
  labs( title = "c)" , x= "Standard Deviation of Lidar Return Intensity", y= "Mean %C")+
  theme_minimal(base_size = 24)
p8<-plot(predict_response(best_model, terms= "zkurt_3")) +
  labs( title = "d)" , x= "Kurtosis of Canopy Height Distribution", y= "Mean %C")+
  theme_minimal(base_size = 24)
p9<-plot(predict_response(best_model, terms= "p1th_3")) +
  labs( title = "e)" , x= "Percentage 1st Returns", y= "Mean %C")+
  theme_minimal(base_size = 24)

p<- (p5+p6 )/(p7 +p8)/( p9)

p
ggsave('plots/marginal_effects.png', p, width = 20, height = 18, dpi = 300)

