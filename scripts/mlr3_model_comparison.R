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
install.packages("genalg")
future::plan("multisession", workers = 5)

#### set up tasks ####
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

###sort out incomplete rows

##Create a task for each data type level  
tsk_all_data = as_task_regr_st(converted_data, target = "wmean_percC",
                               id = "soils_all")


### encoding ####

##Create a task for each data type level  
tsk_all_data = as_task_regr_st(converted_data, target = "wmean_percC",
                               id = "soils_all")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_all_data))[[1]]

#### feature selection ####
instance = fsi(
  task = encoded_tsk,
  learner = lrn("regr.xgboost"),
  resampling = rsmp("spcv_coords", folds = 32),
  measures = msr("regr.rmse"), ## multiple measures works for multicrit but try tpo evaluate different things 
  terminator = trm("evals", n_evals = 100)
)

# Choose optimization algorithm
fselector = fs("genetic_search")

# Run feature selection
progressr::with_progress(fselector$optimize(instance))

all_vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)
  
#### HPO #####
tuning_instance <- ti(
  task = task_st_filter,
  learner = lrn("regr.xgboost"),
  resampling = rsmp("spcv_coords", folds = 30),
  measure = msr("regr.rmse"),
  search_space = lts("regr.xgboost.default"),  # Corrected argument name
  terminator = trm("evals", n_evals = 100)  # Define termination criterion
)

tuner = tnr("mbo")

progressr::with_progress(tuner$optimize(tuning_instance))

#### tuned learner ####
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)
task_st_filter_all <- task_st_filter$clone(deep = TRUE)

#### resampling instance and plot  #####
resampling_instance = mlr3::resample(
  task= task_st_filter_all,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))

dot_plot2 <- function(x, labs = NULL, dot_alpha = 0.01, ...) {
  rr <- x$aggregate(measure = c(
    mlr3::msr("regr.bias"),
    mlr3::msr("regr.rmse"),
    mlr3::msr("regr.mse")
  )) |>
    round(2)
  
  df <- x$prediction() |>
    data.table::as.data.table()
  
  
  # get max of x and y
  max_xy <- round(max(
    max(df$response, na.rm = TRUE),
    max(df$truth, na.rm = TRUE)
  ) / 10) * 10
  
  # plot the results.
  p <- df |>
    ggplot() +
    aes(y = truth, x = response) +
    geom_point(col = "#f08a46", alpha = dot_alpha) +
    #geom_density_2d(aes(col = after_stat(level))) +
    scale_color_viridis_c(direction = -1, option = "mako") +
    guides(alpha = "none", color = "none") +
    geom_abline(slope = 1) +
    #coord_fixed(xlim = c(0, max_xy), ylim = c(0, max_xy)) +
    annotate("text",
             x = max_xy * 0.1, y = max_xy * 0.9,
             label = paste0("bias = ", rr["regr.bias"])
    ) +
    annotate("text",
             x = max_xy * 0.1, y = max_xy * 0.85,
             label = paste0("rmse = ", rr["regr.rmse"])
    ) +
    annotate("text",
             x = max_xy * 0.1, y = max_xy * 0.8,
             label = paste0("mse = ", rr["regr.mse"])
    ) +
    theme_linedraw()
  
  if (!is.null(labs)) {
    p <- p + labs
  }
  
  return(p)
}

P_all<-dot_plot2(resampling_instance, dot_alpha = 0.9)


##### repeat for RS data #####

##Create a task for each data type level  
RS_data_subset<- converted_data%>% select(ends_with("_RS"), (ends_with("_FA")), wmean_percC)
tsk_RS_data = as_task_regr_st(RS_data_subset, target = "wmean_percC",
                              id = "soils_RS")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_RS_data))[[1]]

#### feature selection ####
instance = fsi(
  task = encoded_tsk,
  learner = lrn("regr.xgboost"),
  resampling = rsmp("spcv_coords", folds = 32),
  measures = msr("regr.rmse"), ## multiple measures works for multicrit but try tpo evaluate different things 
  terminator = trm("evals", n_evals = 100)
)

# Choose optimization algorithm
fselector = fs("genetic_search")

# Run feature selection
progressr::with_progress(fselector$optimize(instance))

RS_Vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

#### HPO #####
tuning_instance <- ti(
  task = task_st_filter,
  learner = lrn("regr.xgboost"),
  resampling = rsmp("spcv_coords", folds = 30),
  measure = msr("regr.rmse"),
  search_space = lts("regr.xgboost.default"),  # Corrected argument name
  terminator = trm("evals", n_evals = 100)  # Define termination criterion
)

tuner = tnr("mbo")

progressr::with_progress(tuner$optimize(tuning_instance))

#### tuned learner ####
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_RS <- task_st_filter$clone(deep = TRUE)
#### resampling instance and plot  #####
resampling_instance = mlr3::resample(
  task= task_st_filter_RS,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_RS<-dot_plot2(resampling_instance, dot_alpha = 0.9)


##### repeat for FA data #####

##Create a task for each data type level  
FA_data_subset<- converted_data%>% select((ends_with("_FA")), wmean_percC)
tsk_FA_data = as_task_regr_st(FA_data_subset, target = "wmean_percC",
                              id = "soils_FA")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_FA_data))[[1]]

#### feature selection ####
instance = fsi(
  task = encoded_tsk,
  learner = lrn("regr.xgboost"),
  resampling = rsmp("spcv_coords", folds = 32),
  measures = msr("regr.rmse"), ## multiple measures works for multicrit but try tpo evaluate different things 
  terminator = trm("evals", n_evals = 100)
)

# Choose optimization algorithm
fselector = fs("genetic_search")

# Run feature selection
progressr::with_progress(fselector$optimize(instance))

FA_vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

#### HPO #####
tuning_instance <- ti(
  task = task_st_filter,
  learner = lrn("regr.xgboost"),
  resampling = rsmp("spcv_coords", folds = 30),
  measure = msr("regr.rmse"),
  search_space = lts("regr.xgboost.default"),  # Corrected argument name
  terminator = trm("evals", n_evals = 100)  # Define termination criterion
)

tuner = tnr("mbo")

progressr::with_progress(tuner$optimize(tuning_instance))

#### tuned learner ####
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_FA <- task_st_filter$clone(deep = TRUE)
#### resampling instance and plot  #####
resampling_instance = mlr3::resample(
  task= task_st_filter_FA,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_FA<-dot_plot2(resampling_instance, dot_alpha = 0.9)
#### play around with benchmarking #####

P_all
P_RS
P_FA

tasks = c(task_st_filter_FA, task_st_filter_RS, task_st_filter_all)
learners = lrns(c("regr.xgboost", "regr.featureless"), predict_type = "response")
rsmp_ho = rsmp("spcv_coords", folds = 30)

design = benchmark_grid(tasks, learners, rsmp_ho)
head(design)
bmr = benchmark(design)
bmr$score()

summary_table <- table %>%
  group_by(nr) %>%
  summarise(
    across(where(is.numeric), mean, .names = "wmean_{.col}"),
    across(where(is.character), get_mode, .names = "mode_{.col}"),
    .groups = "drop"
  )

