#### practice #####
install.packages("mlr3spatiotempcv")
install.packages("xgboost")
library(xgboost)
future::plan(future::multisession, workers = future::availableCores()-2)
library(mlr3spatiotempcv)
set.seed(42)
samples_metrics<- read_sf("Data/soil_samples_w_complete_metrics.fgb")
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
  mutate(weight = ifelse(Profundidade_cm_ == "10-20"| Profundidade_cm_ == "20-30", 2, 1))

# Summarize table
summary_df <- samples_metrics %>%
  group_by(Codigo) %>%
  summarise(
    across(where(is.numeric), ~ weighted_mean_func(., weight), .names = "wmean_{.col}"),
    across(where(is.character) & !all_of("Profundidade_cm_"), get_mode, .names = "mode_{.col}"),
    .groups = "drop"
  )

# Remove weight column from summary
summary_df <- summary_df %>%
  select(-contains("weight"))

##drop useless columns
non_sf<-summary_df%>%
  #st_drop_geometry()%>%
  select(!dplyr::any_of(c("Codigo", "wmean_n_RS", "wmean_GEDI.footprints","wmean_id_clean","mode_date_time_FA" )))

## make everything numeric or factor ##
convert_to_numeric_or_factor <- function(df) {
  df %>%
    mutate(across(1:last_col(1), ~ {
      if (suppressWarnings(all(!is.na(as.numeric(.))))) {
        as.numeric(.)
      } else {
        as.factor(.)
      }
    }))
}
converted_data <- convert_to_numeric_or_factor(non_sf)
##make sure headers work in mlr3 ecosystem
#converted_data<- non_sf
colnames(converted_data)<-make.names(colnames(converted_data))

###sort out incomplete rows

##Create a task for each data type level  
tsk_all_data = as_task_regr_st(converted_data, target = "wmean_percC",
                               id = "soils_all")
RS_data_subset<- converted_data%>% select(ends_with("_RS"), (ends_with("_FA")), wmean_percC)
tsk_RS_data = as_task_regr_st(RS_data_subset, target = "wmean_percC",
                            id = "soils_RS")
FA_data_subset<- converted_data%>% select((ends_with("_FA")), wmean_percC)
tsk_FA_data = as_task_regr_st(FA_data_subset, target = "wmean_percC",
                           id = "soils_FA")
#tsk_FA_data$set_col_roles("Profundidade_cm_", roles = "group")

## feature selection 
#TODO create a feature selection sandwich 
# Define filters to use
# Define function to select top-ranked features
select_top_features <- function(task, n_top = 30) {
  
  # Define filters to use
  filters <- c("jmim", "mrmr", "cmim")
  
  # Apply each filter and store results
  filter_results <- lapply(filters, function(flt_name) {
    filter <- flt(flt_name)   # Initialize filter
    filter$calculate(task)    # Compute scores
    data.frame(feature = names(filter$scores), 
               score = rescale(filter$scores), # Normalize scores (0-1)
               filter = flt_name)
  })
  
  # Combine results into one table
  all_scores <- bind_rows(filter_results)
  
  # Aggregate scores (mean score across filters)
  final_ranking <- all_scores %>%
    group_by(feature) %>%
    summarise(agg_score = mean(score)) %>%
    arrange(desc(agg_score))  # Rank from highest to lowest
  
  # Select top `n_top` ranked features
  top_features <- final_ranking %>%
    slice_head(n = min(n_top, nrow(final_ranking))) %>%  # Avoid exceeding available features
    pull(feature)  # Extract feature names
  
  # Keep only the selected features + target variable
  filtered_data <- task$select(top_features)
  
  # Return list with selected features and filtered dataset
  return(list(
    selected_features = top_features,
    filtered_data = filtered_data,
    scores= final_ranking
  ))
}



# Apply function to select top 3 features (for demonstration)
all_result <- select_top_features(tsk_all_data, n_top = 30)
RS_result<-select_top_features(tsk_RS_data, n_top = 30)
FA_result<- select_top_features(tsk_FA_data, n_top = 30)
# Print selected features
print(FA_result$scores)


#TODO hyper parameter tuning spaces hughs or defaults or chat 
#TODO set up work in pipeline 
lts_rpart = lts("regr.xgboost.default")
instance = ti(
  task = tsk_FA_data,
  learner = lrn("regr.xgboost"),
  resampling = rsmp("cv", folds = 3),
  measures = msr("regr.mae"),
  terminator = trm("evals", n_evals = 20),
  search_space = lts_rpart
)

RF_learner<-lrn("regr.ranger")
RS_instance = fselect(
  fselector = fs("sequential"),
  task =  tsk_RS_data,
  learner = RF_learner,
  resampling = rsmp("holdout"),
  measure = msr("regr.mae")
)
all_instance = fselect(
  fselector = fs("sequential"),
  task =  tsk_all_data,
  learner = RF_learner,
  resampling = rsmp("holdout"),
  measure = msr("regr.mae")
)
FA_instance = fselect(
  fselector = fs("sequential"),
  task =  tsk_FA_data,
  learner = RF_learner,
  resampling = rsmp("holdout"),
  measure = msr("regr.mae")
)
dt = as.data.table(instance$archive)

autoplot(FA_instance, type = "performance")
autoplot(all_instance, type = "performance")
autoplot(RS_instance, type = "performance")

## set up models based on feature selection ( ideal combo from sequential plus top 3 IG)
RS_keep = paste0(c(RS_instance$result_feature_set, names(head(RS_flt_gain$scores, 3))))
all_keep = paste0(c(all_instance$result_feature_set, names(head(all_flt_gain$scores, 3))))
FA_keep = paste0(c(FA_instance$result_feature_set, names(head(FA_flt_gain$scores, 3))))
tsk_all_data$select(all_keep)

tsk_RS_data$select(RS_keep)

tsk_FA_data$select(FA_keep)

## use benchmark to compare learners
tasks = c(FA_result, all_result, RS_result)
learners = lrns(c("regr.ranger", "regr.featureless"), predict_type = "response")
rsmp_ho = rsmp("holdout")

design = benchmark_grid(tasks, learners, rsmp_ho)
head(design)
bmr = benchmark(design)
bmr$score()



###### figure out model tuning #######
# Load required libraries
library(mlr3)
library(mlr3learners)
library(mlr3filters)
library(mlr3tuning)
library(mlr3pipelines)
library(mlr3mbo)   # Bayesian optimization
library(paradox)   # Hyperparameter search space
library(dplyr)
library(scales)
library(mlr3fselect)

##multicrypt

# Define function to select top-ranked features
select_top_features <- function(task, n_top = 30) {
  
  # Define filters to use
  filters <- c("jmim", "mrmr", "cmim")
  
  # Apply each filter and store results
  filter_results <- lapply(filters, function(flt_name) {
    filter <- flt(flt_name)   # Initialize filter
    filter$calculate(task)    # Compute scores
    data.frame(feature = names(filter$scores), 
               score = rescale(filter$scores), # Normalize scores (0-1)
               filter = flt_name)
  })
  
  # Combine results into one table
  all_scores <- bind_rows(filter_results)
  
  # Aggregate scores (mean score across filters)
  final_ranking <- all_scores %>%
    group_by(feature) %>%
    summarise(agg_score = mean(score)) %>%
    arrange(desc(agg_score))  # Rank from highest to lowest
  
  # Select top `n_top` ranked features
  top_features <- final_ranking %>%
    slice_head(n = min(n_top, nrow(final_ranking))) %>%
    pull(feature)  # Extract feature names
  
  # Keep only the selected features + target variable
  filtered_data <- task$data(cols = c(top_features, task$target_names))
  
  # Return list with selected features and filtered dataset
  return(list(
    selected_features = top_features,
    filtered_data = filtered_data
  ))
}

# Load your regression task (replace `tsk_FA_data` with actual task)
tsk_all_data = as_task_regr_st(converted_data, target = "wmean_percC",
                              id = "soils_all")  # Set target variable

# Perform initial feature selection
feature_selection_result <- select_top_features(tsk_all_data, n_top = 30)

# Create a new task with selected features
tsk_filtered <- TaskRegr$new(id = "filtered_task", 
                             backend = feature_selection_result$filtered_data, 
                             target = "wmean_percC")

# Define learners
learners <- list(
  rpart = lrn("regr.rpart"),
  ranger = lrn("regr.ranger")
)

# Define hyperparameter search spaces
param_sets <- list(
  rpart = lts("regr.rpart.default"),
  
  ranger = lts("regr.ranger.default")
)

# Define Bayesian Optimizer
bayes_optim <- function(learner, param_set, task) {
  instance <- TuningInstanceBatchSingleCrit$new(
    task = task,
    learner = learner,
    resampling = rsmp("holdout"),
    measure = msr("regr.rmse"),
    search_space = param_set,  # Corrected argument name
    terminator = trm("evals", n_evals = 20)  # Define termination criterion
  )
  
  tuner <- TunerMbo$new()  # Initialize Bayesian Optimizer
  tuner$optimize(instance) # Run optimization
  
  return(instance$result$params)
}

# Run optimization for both learners
opt_results <- list()
for (name in names(learners)) {
  opt_results[[name]] <- bayes_optim(learners[[name]], param_sets[[name]], tsk_filtered)
}

# Train final models with optimized hyperparameters
trained_models <- list()
for (name in names(learners)) {
  learner <- learners[[name]]
  
  # Set optimized parameters
  learner$param_set$values <- opt_results[[name]]
  

  
  # Perform cross-validation (3-fold)
  resampling <- rsmp("cv", folds = 3)
  resample_result <- mlr3::resample(tsk_filtered, learner, resampling)
  
  # Store results
  trained_models[[name]] <- list(
    model = learner,
    performance = resample_result$aggregate(msrs(c("regr.rmse", "regr.mae")))
  )
}

# Compare model performance
performance_comparison <- data.frame(
  Model = c("rpart", "Ranger"),
  RMSE = c(trained_models$rpart$performance[1], 
           trained_models$ranger$performance[1]),
  MAE = c(trained_models$rpart$performance[2], 
           trained_models$ranger$performance[2])
)

# Print final comparison
print(performance_comparison)

 

##### apply settings to model #####
# Retrieve best hyperparameters for ranger
best_params <- opt_results[["ranger"]]

# Create a new ranger learner and set the best parameters
final_learner <- lrn("regr.ranger")
final_learner$param_set$values <- best_params
# Enable feature importance calculation
final_learner$param_set$values$importance <- "permutation"
# Set up cross-validation
resampling <- rsmp("cv", folds = 3)

# Run cross-validation
final_resample <- mlr3::resample(tsk_filtered, final_learner, resampling)
# Aggregate results across CV folds
final_performance <- final_resample$aggregate(msrs(c("regr.rmse", "regr.mae")))

# Print the final model's performance
print(final_performance)
# Train the final model on the full dataset
final_learner$train(tsk_filtered)

# Extract and display feature importance
importance <- final_learner$importance()
print(importance)
# Get predictions on the full dataset
predictions <- final_learner$predict(tsk_filtered)

# Convert to data frame
pred_df <- data.frame(
  Actual = predictions$truth,
  Predicted = predictions$response
)

# Scatter plot: Actual vs. Predicted
ggplot(pred_df, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + # Ideal line
  theme_minimal() +
  labs(title = "Predicted vs. Actual Values",
       x = "Actual Values", y = "Predicted Values")
