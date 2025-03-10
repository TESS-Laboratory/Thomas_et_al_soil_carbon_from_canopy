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

set.seed(123)  # For reproducibility

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


#### Function to Split Dataset (Preserving sf Structure) ####
split_dataset <- function(data, split_ratio = 0.7) {
  if (!inherits(data, "sf")) {
    stop("Error: Input data must be an sf object!")
  }
  
  total_rows <- nrow(data)
  if (total_rows < 2) {
    stop("Error: Not enough data for splitting!")
  }
  
  train_indices <- sample(seq_len(total_rows), size = floor(split_ratio * total_rows))
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  return(list(train = train_data, test = test_data))
}

### Function to Run the Spatial Machine Learning Pipeline ####
run_ml_pipeline <- function(train_data, test_data, feature_suffixes, folds_spcv = 5, folds_hpo = 5, sim_id) {
  message(paste0("Running Simulation ", sim_id, " with features: ", paste(feature_suffixes, collapse = ", ")))
  
  # Extract Selected Features
  #selected_features <- select("wmean_percC_5", ends_with(feature_suffixes))  
  train_data_filtered <- train_data %>% select("wmean_percC_5", ends_with(feature_suffixes))
  test_data_filtered <- test_data %>% select("wmean_percC_5", ends_with(feature_suffixes))
  
  # Convert to Regression Task (for Spatial Data)
  task <- as_task_regr_st(train_data_filtered, target = "wmean_percC_5", id = paste0("sim_", sim_id))
  
  poe = po("encode")
  encoded_tsk = poe$train(list(task))[[1]]
  
  ## feature selection ##
  instance = fsi(
    task = encoded_tsk,
    learner = lrn("regr.xgboost"),
    resampling = rsmp("spcv_coords", folds = 10),
    measures = msr("regr.rmse"), ## multiple measures works for multicrit but try tpo evaluate different things 
    terminator = trm("evals", n_evals = 10)
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
    learner = lrn("regr.xgboost"),
    resampling = rsmp("spcv_coords", folds = 10),
    measure = msr("regr.rmse"),
    search_space = lts("regr.xgboost.default"),  # Corrected argument name
    terminator = trm("evals", n_evals = 10)  # Define termination criterion
  )
  
  tuner = tnr("mbo")
  
  progressr::with_progress(tuner$optimize(tuning_instance))
  
  ## tuned learner ##
  tuned_learner = lrn("regr.xgboost")
  tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)
  
  task_st_filter_S1 <- task_st_filter$clone(deep = TRUE)
  ## resampling instance and plot  ##
  resampling_instance = mlr3::resample(
    task= task_st_filter_S1,
    learner = tuned_learner,
    resampling = rsmp("spcv_coords", folds = 10))
  
  learner_final <- tuned_learner$train(task_st_filter_S1)
  predictions <- learner_final$predict_newdata(test_data_filtered)
  
  # 5. Compute Metrics
  rmse <- sqrt(mean((predictions$response - test_data_filtered$wmean_percC_5)^2))
  mae <- mean(abs(predictions$response - test_data_filtered$wmean_percC_5))
  bias <- mean(predictions$response - test_data_filtered$wmean_percC_5)
  
  metrics <- data.frame(sim_id = sim_id, RMSE = rmse, MAE = mae, Bias = bias)
  
  # 6. Modelled vs Observed Plot
  plot <- ggplot(test_data_filtered, aes(x = wmean_percC_5, y = predictions$response)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = paste0("Simulation", sim_id),
      x = "Observed Carbon",
      y = "Modelled Carbon"
    ) +
    annotate("text", x = min(test_data_filtered$wmean_percC_5), y = max(predictions$response)*0.97 , 
             label = paste0("RMSE:", round(rmse, 2),"\n",
                            "MAE:", round(mae, 2),"\n",
                            "Bias:", round(bias, 2)),
             hjust = 0.1, size = 2, color = "black") +
    theme_minimal()
  
  return(list(rr = resampling_instance, metrics = metrics, plot = plot, retained_variables = S1_Vars))
}


# Split Data into Train and Test
split_data <- split_dataset(converted_data)
train_data <- split_data$train
test_data <- split_data$test

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
  list(feature_suffixes = c("_1", "_2", "_3"), sim_id = 10)
)

# Run Simulations
results <- lapply(simulations, function(sim) {
  run_ml_pipeline(train_data, test_data, sim$feature_suffixes, sim_id = sim$sim_id)
})

# Collect Benchmark Metrics
metrics_df <- do.call(rbind, lapply(results, function(res) res$metrics))

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
