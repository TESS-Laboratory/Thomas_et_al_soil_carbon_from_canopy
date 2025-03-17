## fully modular mlr3 script
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
install.packages("remotes")
remotes::install_github("mlr-org/mlr3extralearners@*release")
install.packages("glmnet")
install.packages("kknn")
library(glmnet)
library(kknn)
library(mlr3extralearners)
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

### read in data and set parameters ####

#dataset
samples_metrics<- read_sf("Data/soil_samples_w_complete_metrics.fgb")

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

# Define learners and search space Configurations
students <- list(
  list(learner = "regr.ranger", SS = "regr.ranger.default"),
  #list(learner = "regr.glm", SS = "regr.glmnet.default"),
  #list(learner = "regr.kknn", SS = "regr.kknn.default"),
  #list(learner = "regr.km", SS = "regr.ranger.default"),
  #list(learner = "regr.nnet", SS = "regr.ranger.default"),
  list(learner = "regr.rpart", SS = "regr.rpart.default"),
  #list(learner = "regr.svm", SS = "regr.svm.default"),
  list(learner = "regr.xgboost", SS = "regr.xgboost.rbv1")
)


### create task list ####
set_up_tasks <- function(train_data, feature_suffixes, sim_id) {
  message(paste0("Running Simulation ", sim_id, " with features: ", paste(feature_suffixes, collapse = ", ")))
  
  # Extract Selected Features
  train_data_filtered <- train_data %>% select("wmean_percC_5", ends_with(feature_suffixes))
  
  # Convert to Regression Task (for Spatial Data)
  task <- as_task_regr_st(train_data_filtered, target = "wmean_percC_5", id = paste0("sim_", sim_id))
  poe = po("encode")
  encoded_tsk = poe$train(list(task))[[1]]
  
  return(encoded_tsk)  
}

task_list <- lapply(simulations, function(sim) {
  set_up_tasks(converted_data, sim$feature_suffixes, sim_id = sim$sim_id)
})

### create at function for multiple at learners ####
set_up_at_learners <- function(lrnr, SS) {
  
  at = auto_tuner(
    tuner = tnr("mbo"),
    learner = lrn(lrnr),
    resampling = rsmp("spcv_coords", folds = 3),
    measure = msr("regr.rmse"),
    search_space = lts(SS),
    term_evals = 10,
    terminator = trm("evals", n_evals = 100)
  )
  
  afs = auto_fselector(
    fselector = fs("genetic_search"),
    learner = at,
    resampling = rsmp("spcv_coords", folds = 3),
    measure = msr("regr.rmse"),
    term_evals = 10,
    terminator = trm("evals", n_evals = 100)
  )
  
  return( afs)  
}

learns_list <- lapply(students, function(student) {
  set_up_at_learners(lrnr = student$learner, SS = student$SS)
})
### bench mark design ####

tasks = task_list
learners = learns_list
resamplings = rsmp("spcv_coords", folds = 3)

design = benchmark_grid(task_list, learners, resamplings)
bmr = benchmark(design)
bmr$aggregate

### extract results ####
best_row <- bmr$aggregate()[which.min(bmr$aggregate()$regr.mse), ]
learn <- best_row$learner_id

result <- bmr$aggregate %>%
  group_by(task_id) %>%
  filter(reg.mse == min(reg.mse, na.rm = TRUE)) %>%
  ungroup()

### resample best learners ####


### plots and analysis ####