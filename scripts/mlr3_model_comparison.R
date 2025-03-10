# Load required libraries
library(mlr3)
library(patchwork)
library(mlr3spatiotempcv)
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
library(genalg)
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
  select(!dplyr::any_of(c("Codigo", "wmean_n_3", "wmean_GEDI.footprints","wmean_id_clean","mode_date_time_FA" )))

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





### encoding ####

##Create a task for each data type level  
tsk_S1 = as_task_regr_st(converted_data, target = "wmean_percC_5",
                         id = "sim_1")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S1))[[1]]

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

##### repeat for sim 1 #####


S1_data_subset<- converted_data%>% select(ends_with("_1"),  wmean_percC_5)
tsk_S1 = as_task_regr_st(S1_data_subset, target = "wmean_percC_5",
                         id = "sim_1")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S1))[[1]]

## feature selection ##
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
  resampling = rsmp("spcv_coords", folds = 30),
  measure = msr("regr.rmse"),
  search_space = lts("regr.xgboost.default"),  # Corrected argument name
  terminator = trm("evals", n_evals = 100)  # Define termination criterion
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
  resampling = rsmp("spcv_coords", folds = 30))


P_S1<-dot_plot2(resampling_instance, dot_alpha = 0.9)


##### repeat for sim 2 #####

 
S2_data_subset<- converted_data%>% select(ends_with(c("_1", "_2")), wmean_percC_5)
tsk_S2 = as_task_regr_st(S2_data_subset, target = "wmean_percC_5",
                              id = "sim_2")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S2))[[1]]

## feature selection ##
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

S2_Vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

## HPO ##
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

## tuned learner ##
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_S2 <- task_st_filter$clone(deep = TRUE)
## resampling instance and plot  ##
resampling_instance = mlr3::resample(
  task= task_st_filter_S2,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_S2<-dot_plot2(resampling_instance, dot_alpha = 0.9)

##### repeat for sim 3 #####


S3_data_subset<- converted_data%>% select(ends_with(c("_3")), wmean_percC_5)
tsk_S3 = as_task_regr_st(S3_data_subset, target = "wmean_percC_5",
                         id = "sim_3")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S3))[[1]]

## feature selection ##
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

S3_Vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

## HPO ##
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

## tuned learner ##
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_S3 <- task_st_filter$clone(deep = TRUE)
## resampling instance and plot  ##
resampling_instance = mlr3::resample(
  task= task_st_filter_S3,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_S3<-dot_plot2(resampling_instance, dot_alpha = 0.9)

##### repeat for sim 4 #####


S4_data_subset<- converted_data%>% select(ends_with(c("_5")), wmean_percC_5)
tsk_S4 = as_task_regr_st(S4_data_subset, target = "wmean_percC_5",
                         id = "sim_4")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S4))[[1]]

## feature selection ##
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

S4_Vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

## HPO ##
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

## tuned learner ##
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_S4 <- task_st_filter$clone(deep = TRUE)
## resampling instance and plot  ##
resampling_instance = mlr3::resample(
  task= task_st_filter_S4,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_S4<-dot_plot2(resampling_instance, dot_alpha = 0.9)

##### repeat for sim 5 #####


S5_data_subset<- converted_data%>% select(ends_with(c("_4", "_5")), wmean_percC_5)
tsk_S5 = as_task_regr_st(S5_data_subset, target = "wmean_percC_5",
                         id = "sim_5")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S5))[[1]]

## feature selection ##
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

S5_Vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

## HPO ##
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

## tuned learner ##
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_S5 <- task_st_filter$clone(deep = TRUE)
## resampling instance and plot  ##
resampling_instance = mlr3::resample(
  task= task_st_filter_S5,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_S5<-dot_plot2(resampling_instance, dot_alpha = 0.9)

##### repeat for sim 6 #####


S6_data_subset<- converted_data%>% select(ends_with(c("_1", "_5")), wmean_percC_5)
tsk_S6 = as_task_regr_st(S6_data_subset, target = "wmean_percC_5",
                         id = "sim_6")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S6))[[1]]

## feature selection ##
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

S6_Vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

## HPO ##
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

## tuned learner ##
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_S6 <- task_st_filter$clone(deep = TRUE)
## resampling instance and plot  ##
resampling_instance = mlr3::resample(
  task= task_st_filter_S6,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_S6<-dot_plot2(resampling_instance, dot_alpha = 0.9)

##### repeat for sim 7 #####


S7_data_subset<- converted_data%>% select(ends_with(c("_1", "_2", "_4", "_5")), wmean_percC_5)
tsk_S7 = as_task_regr_st(S7_data_subset, target = "wmean_percC_5",
                         id = "sim_7")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S7))[[1]]

##feature selection ##
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

S7_Vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

## HPO ##
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

## tuned learner ##
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_S7 <- task_st_filter$clone(deep = TRUE)
## resampling instance and plot  ##
resampling_instance = mlr3::resample(
  task= task_st_filter_S7,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_S7<-dot_plot2(resampling_instance, dot_alpha = 0.9)

##### repeat for sim 8 #####


S8_data_subset<- converted_data%>% select(ends_with(c("_3", "_4", "_5")), wmean_percC_5)
tsk_S8 = as_task_regr_st(S8_data_subset, target = "wmean_percC_5",
                         id = "sim_8")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S8))[[1]]

## feature selection ##
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

S8_Vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

## HPO ##
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

## tuned learner ##
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_S8 <- task_st_filter$clone(deep = TRUE)
## resampling instance and plot  ##
resampling_instance = mlr3::resample(
  task= task_st_filter_S8,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_S8<-dot_plot2(resampling_instance, dot_alpha = 0.9)

##### repeat for sim 9 #####


tsk_S9 = as_task_regr_st(converted_data, target = "wmean_percC_5",
                         id = "sim_9")
poe = po("encode")
encoded_tsk = poe$train(list(tsk_S9))[[1]]

# feature selection #
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

S9_Vars<-instance$result$features
## create new feature
task_st_filter <- encoded_tsk$clone(deep = TRUE)
task_st_filter$select(
  instance$result$features[[1]]
)

# HPO #
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

# tuned learner #
tuned_learner = lrn("regr.xgboost")
tuned_learner$param_set$set_values(.values= tuning_instance$result_learner_param_vals)

task_st_filter_S9 <- task_st_filter$clone(deep = TRUE)
# resampling instance and plot  #
resampling_instance = mlr3::resample(
  task= task_st_filter_S9,
  learner = tuned_learner,
  resampling = rsmp("spcv_coords", folds = 30))


P_S9<-dot_plot2(resampling_instance, dot_alpha = 0.9)
