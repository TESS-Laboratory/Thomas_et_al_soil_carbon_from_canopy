### set up environment ####
install.packages("mlr3")
install.packages("mlr3spatiotempcv")
install.packages("mlr3learners")
install.packages("mlr3filters")
install.packages("mlr3fselect")
install.packages("mlr3tuning")
install.packages("mlr3mbo")
install.packages("data.table")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("genalg")
install.packages("xgboost")
install.packages("DiceKriging")
install.packages("rgenoud")
install.packages("ranger")
install.packages("mlr3pipelines")
install.packages("paradox")
install.packages("scales")
install.packages("progressr")
install.packages("mlr3verse")
install.packages("mlr3tuningspaces")
install.packages("FSelectorRcpp")
library(FSelectorRcpp)
library(mlr3tuningspaces)
library(mlr3verse)
library(progressr)
library(scales)
library(paradox)
library(mlr3pipelines)
library(ranger)
library(rgenoud)
library(DiceKriging)
library(xgboost)
library(genalg)
library(mlr3)
library(mlr3spatiotempcv)
library(mlr3learners)
library(mlr3filters)
library(mlr3fselect)
library(mlr3tuning)
library(mlr3mbo)
library(data.table)
library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)  # For arranging plots





set.seed(42)  # For reproducibility

samples_metrics<- read_sf("Data/soil_samples_w_complete_metrics.fgb")

#### clean #####

## take averages for soil depth
weighted_mean_func <- function(x, w) {
  if (is.numeric(x)) {
    return(weighted.mean(x, w, na.rm = TRUE))
  } else {
    return(NA)  # Placeholder for non-numeric columns
  }
}

# Define a function to calculate mode
get_mode <- function(x) {
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]
}

# Compute weights: "yes" = 2, "no" = 1
samples_metrics <- samples_metrics %>%
  mutate(weight = ifelse(Profundidade_cm__5 == "10-20"| Profundidade_cm__5 == "20-30", 2, 1))

# Summarize table
summary_df <- samples_metrics %>%
  group_by(Codigo) %>%
  summarise(
    across(where(is.numeric), ~ weighted_mean_func(., weight), .names = "wmean_{.col}"),
    across(where(is.character) & !all_of("Profundidade_cm__5"), get_mode, .names = "mode_{.col}"),
    .groups = "drop"
  )

# Remove weight column from summary
summary_df <- summary_df %>%
  select(-contains(c("weight", "CperN")))

##drop useless columns
non_sf<-summary_df%>%
  #st_drop_geometry()%>%
  select(!dplyr::any_of(c("Codigo", "wmean_n_3", "wmean_GEDI.footprints","wmean_id_clean","mode_date_time_2" )))

## make everything numeric or factor ##
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

converted_data <- convert_to_numeric_or_factor(non_sf)
##make sure headers work in mlr3 ecosystem
#converted_data<- non_sf
colnames(converted_data)<-make.names(colnames(converted_data))

### start modelling ####


# Function to Run the Spatial Machine Learning Pipeline
run_ml_pipeline <- function(train_data, feature_suffixes, folds_spcv = 6, folds_hpo = 6, sim_id, xy_lim= 8) {
  message(paste0("Running Simulation ", sim_id, " with features: ", paste(feature_suffixes, collapse = ", ")))
  
  # Extract Selected Features
  train_data_filtered <- train_data %>% select("wmean_percC_5", ends_with(feature_suffixes))
  
  
  
  
  # Convert to Regression Task (for Spatial Data)
  task <- as_task_regr_st(train_data_filtered, target = "wmean_percC_5", id = paste0("sim_", sim_id))
  
  # select best learner
  learners = lrns(c( "regr.ranger", "regr.rpart", "regr.xgboost"), predict_type = "response")
  rsmp_cv5 = rsmp("cv", folds = 5)
  design = benchmark_grid(task, learners, rsmp_cv5)
  bmr = benchmark(design)
  
  best_row <- bmr$aggregate()[which.min(bmr$aggregate()$regr.mse), ]
  learn <- best_row$learner_id
  
  poe = po("encode")
  encoded_tsk = poe$train(list(task))[[1]]
  
  ## feature selection ##
  instance = fsi(
    task = encoded_tsk,
    learner = lrn(learn),
    resampling = rsmp("spcv_coords", folds = folds_spcv),
    measures = msr("regr.rmse"), ## multiple measures works for multicrit but try tpo evaluate different things 
    terminator = trm("evals", n_evals = 100)
  )
  
  # Choose optimization algorithm
  fselector = fs("genetic_search")
  
  # Run feature selection
  progressr::with_progress(fselector$optimize(instance))
  
  S1_Vars<-instance$result$features
  ## create new feature
  task_st_filter <- encoded_tsk$clone(deep = TRUE)
  task_st_filter$select(
    instance$result$features[[1]]
  )
  
  ## HPO ##
  tuning_instance <- ti(
    task = task_st_filter,
    learner = lrn(learn),
    resampling = rsmp("spcv_coords", folds = folds_hpo),
    measure = msr("regr.rmse"),
    search_space = lts(paste0(learn,".default")),  # Corrected argument name
    terminator = trm("evals", n_evals = 100)  # Define termination criterion
  )
  
  tuner = tnr("mbo")
  
  progressr::with_progress(tuner$optimize(tuning_instance))
  
  ## tuned learner ##
  tuned_learner = lrn(learn)
  tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)
  
  task_st_filter_S1 <- task_st_filter$clone(deep = TRUE)
  ## resampling instance and plot  ##
  resampling_instance = mlr3::resample(
    task= task_st_filter_S1,
    learner = tuned_learner,
    resampling = rsmp("spcv_coords", folds = 10))
  
  learner_final <- tuned_learner$train(task_st_filter_S1)
  
  
  # 5. Compute Metrics
  rr <- resampling_instance$aggregate(measure = c(
    mlr3::msr("regr.bias"),
    mlr3::msr("regr.rmse"),
    mlr3::msr("regr.mse")
  )) |>
    round(2)
  
  df <- resampling_instance$prediction() |>
    data.table::as.data.table()
  
  
  # 6. Modelled vs Observed Plot
  # get max of x and y
  max_xy <- round(max(
    max(df$response, na.rm = TRUE),
    max(df$truth, na.rm = TRUE)
  ) / 10) * 10
  
  # plot the results.
  Plot <- df |>
    ggplot() +
    aes(y = response, x = truth) +
    geom_point(col = "#f08a46", alpha = 0.9) +
    #geom_density_2d(aes(col = after_stat(level))) +
    scale_color_viridis_c(direction = -1, option = "mako") +
    guides(alpha = "none", color = "none") +
    geom_abline(slope = 1) +
    coord_fixed(xlim = c(0, xy_lim), ylim = c(0, xy_lim)) +
    labs(
      title = paste0("Simulation", sim_id),
      x = "Observed Carbon",
      y = "Modelled Carbon"
    ) +
    annotate("text",
             x = xy_lim * 0.15, y = xy_lim * 0.95,
             label = paste0("bias = ", rr["regr.bias"])
    ) +
    annotate("text",
             x = xy_lim * 0.15, y = xy_lim * 0.9,
             label = paste0("rmse = ", rr["regr.rmse"])
    ) +
    annotate("text",
             x = xy_lim * 0.15, y = xy_lim * 0.85,
             label = paste0("mse = ", rr["regr.mse"])
    ) +
    theme_linedraw()
  
  
  
  return(list(rr = resampling_instance, metrics = rr, plot = Plot, retained_variables = S1_Vars, learner = learner_final)) ## add learner to this list 
}




# Define Simulation Configurations
simulations <- list(
  list(feature_suffixes = c("_1"), sim_id = 1),
  list(feature_suffixes = c("_1", "_2"), sim_id = 2),
  list(feature_suffixes = c("_3"), sim_id = 3),
  list(feature_suffixes = c("_5"), sim_id = 4),
  list(feature_suffixes = c("_4", "_5"), sim_id = 5),
  list(feature_suffixes = c("_1", "_5"), sim_id = 6),
  list(feature_suffixes = c("_1", "_2", "_4", "_5"), sim_id = 7),
  list(feature_suffixes = c("_3", "_4", "_5"), sim_id = 8),
  list(feature_suffixes = c("_1", "_2", "_4", "_5", "_3"), sim_id = 9),
  list(feature_suffixes = c("_1", "_2", "_4", "_3"), sim_id = 10),
  list(feature_suffixes = c("_1", "_2", "_3"), sim_id = 11)
)

# Run Simulations
results <- lapply(simulations, function(sim) {
  run_ml_pipeline(converted_data, sim$feature_suffixes, sim_id = sim$sim_id)
})

# Collect Benchmark Metrics
metrics_df <- do.call(rbind, lapply(results, function(res) res$metrics))
variables_df <- do.call(rbind, lapply(results, function(res) res$retained_variables[1]))

# Plot Benchmarking Results
benchmark_plot <- ggplot(metrics_df, aes(x = as.factor(sim_id))) +
  geom_bar(aes(y = RMSE, fill = "RMSE"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = MAE, fill = "MAE"), stat = "identity", position = "dodge") +
  labs(title = "Benchmarking Results", x = "Simulation ID", y = "Error Metrics") +
  scale_fill_manual(values = c("RMSE" = "blue", "MAE" = "red")) +
  theme_minimal()

# Arrange Modelled vs Observed Plots in a 3x3 Grid
plot_list <- lapply(results, function(res) res$plot)
final_grid_plot <- wrap_plots(plot_list, ncol = 3, nrow = 4)

# Display Results
print(benchmark_plot)
print(final_grid_plot)
