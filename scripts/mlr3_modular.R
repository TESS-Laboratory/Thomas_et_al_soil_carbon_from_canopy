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
library(tidyr)
library(ggplot2)
library(sf)
library(patchwork)  # For arranging plots
library(future)

future::plan("multisession", workers = 30)

#' Create hyperparameter tuning space for xgboost
#' @param prefix character prefix to add to parameter names - useful for
#' ensemble models or where the xgboost leanrer id has been changed.
#' @param add_paramsets an additional ParamSet object to add to the
#' collection
#' @return a paradox ParamSet object
#' @export
xgboost_ps <- function(prefix = NULL, add_paramset = NULL) {
  prefix <- if (!is.null(prefix)) paste0(prefix, ".") else ""
  if (!is.null(add_paramset)) paradox::assert_param_set(add_paramset)
  paramset <- setNames(
    list(
      paradox::p_dbl(0.01, 0.5, logscale = TRUE), # eta
      paradox::p_int(1, 1000), # nrounds
      paradox::p_int(3, 12), # max_depth
      paradox::p_dbl(0.1, 1), # subsample
      paradox::p_dbl(0.1, 1), # colsample_bytree
      paradox::p_dbl(0.1, 1), # colsample_bylevel
      paradox::p_dbl(0.001, 700, logscale = TRUE), # alpha
      paradox::p_dbl(0.001, 700, logscale = TRUE) # lambda
    ),
    paste0(
      prefix,
      c(
        "eta",
        "nrounds",
        "max_depth",
        "subsample",
        "colsample_bytree",
        "colsample_bylevel",
        "alpha",
        "lambda"
      )
    )
  )
  
}

set.seed(42)  # For reproducibility
mlr_measures$get("regr.smape")

### read in data and set parameters ####

#dataset
samples_metrics<- read_sf("~/workspace/PhD_work/soil_chapter/Data/soil_samples_w_complete_metrics.fgb")
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

##TODO filter metrics with some kind of correlation metric or most important (20 odd) 

##make sure headers work in mlr3 ecosystem
#converted_data<- non_sf
colnames(converted_data)<-make.names(colnames(converted_data))

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
  list(feature_suffixes = c("_1", "_2", "_3"), sim_id = 11),
  list(feature_suffixes = c("_1", "_3"), sim_id = 12)
)

# Define learners and search space Configurations
students <- list(
  list(learner = lrn("regr.ranger", importance = "impurity"), SS = "regr.ranger.default"),
  list(learner = lrn("regr.glm"), SS=NULL),
  #list(learner = lrn("regr.kknn"), SS = "regr.kknn.default"),
  list(learner = lrn("regr.rpart"), SS = "regr.rpart.default")
  #list(learner = lrn("regr.svm"), SS = "regr.svm.default")
  #list(learner = lrn("regr.xgboost"), SS = xgboost_ps)
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
  if(is.null(SS)){
    at=lrnr
  } else{
    at = auto_tuner(
      tuner = tnr("mbo"),
      learner = lrnr,
      resampling = rsmp("spcv_coords", folds = 3),
      measure = msr("regr.smape"),
      search_space = lts(SS),
      term_evals = 100, ## still might not be enough! 
      terminator = trm("evals", n_evals = 100)
    )
  }
  
  afs = auto_fselector(
    fselector = fs("genetic_search"), ##TODO try forwards selection can't use paralisation but can set up max feature number 
    learner = at,
    resampling = rsmp("spcv_coords", folds = 3),
    measure = msr("regr.smape"),
    term_evals = 100,
    terminator = trm("evals", n_evals = 100)
  )
  
  return(afs)  
}

learns_list <- lapply(students, function(students) {
  set_up_at_learners(lrnr = students$learner, SS = students$SS)
})
### bench mark design ####

tasks = task_list
learners = learns_list
resamplings = rsmp("spcv_coords", folds = 3)

design = benchmark_grid(task_list, learners, resamplings)
bmr = progressr::with_progress(benchmark(design, store_models = TRUE))
bmr$aggregate()

### extract results ####
best_row <- bmr$aggregate()[which.min(bmr$aggregate()$regr.mse), ]
learn <- best_row$learner_id

x <- bmr$aggregate(measures = c(msr("regr.rmse"), msr("regr.mse"))) %>%
  group_by(task_id) %>%
  filter(regr.rmse == min(regr.rmse, na.rm = TRUE)) %>%
  ungroup()

m_obs<- mean(converted_data$wmean_percC_5)
x$rRMSE<- x$regr.rmse/m_obs
xforexport<- x%>% select(-"resample_result")
write.csv(xforexport, "Data/sims_w_measures.csv", row.names = FALSE)
x <- x %>%
  arrange(regr.rmse) %>%
  mutate(task_id = factor(task_id, levels = rev(task_id)))  # reversed for top-down

# Transform to long format (without negating)
df_long <- x %>%
  pivot_longer(cols = c(regr.mse, regr.rmse), names_to = "metric", values_to = "value")

# Plot as grouped horizontal bar chart
p<-ggplot(df_long, aes(x = task_id, y = value, fill = metric)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("regr.mse" = "skyblue", "regr.rmse" = "tomato")) +
  coord_flip() +  # Flip to make it horizontal
  labs(
    y = "Metric Value",
    x = "Model (task_id)",
    title = "Grouped Horizontal Bar Plot of Model Measures"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank()
  )
p
ggsave("Plots/barchart_of_sim_performance.png", plot = p, width = 12, height = 8, dpi = 300)
### resample best learners ####
##bmr$score()$learner[[2]]$importance()

features_function<- function (i){
  features<-extract_inner_fselect_results(bmr)$features[[i]]
  return(features)
}


y<- sort(x$nr)

# y
#[1]  2  5  9 11 14 16 19 22 27 30 31
z<- c(y*3-2,y*3-1, y*3)
#> z
#[1]  4 13 25 31 40 46 55 64 79 88 91  5 14 26 32 41 47 56 65 80 89 92  6 15 27 33 42 48 57 66 81 90 93
retained_features<- lapply(z, features_function)
collapsed_unique <- unique(unlist(retained_features))

# Step 1: Flatten the list
all_values <- unlist(retained_features)

# Step 2: Count frequencies
value_counts <- table(all_values)

# Step 3: Convert to data frame
df <- as.data.frame(value_counts)
colnames(df) <- c("value", "count")

# Step 4: Apply division rules
df$adjusted_count <- with(df, ifelse(
  grepl("_1$", value), count / 8,
  ifelse(grepl("_2$", value), count / 5,
         ifelse(grepl("_3$", value), count / 6,
                ifelse(grepl("_4$", value), count / 5,
                       ifelse(grepl("_5$", value), count / 6, count)
                )))))

df <- df[order(-df$adjusted_count), ]

# Print the sorted table
print(df)
### plots and analysis ####
write.csv(df, "~/workspace/PhD_work/soil_chapter/Data/variable_presence_count.csv", row.names = FALSE)

#autoplot(bmr$score(predictions = TRUE)$prediction_test[[x]])
## grid plot of truth vs response
plot_sims<- function (x, xy_lim= 8) {
  df1 <- bmr$score(predictions = TRUE, ids = TRUE)[x*3]$prediction_test[[1]] |>
    data.table::as.data.table()
  df2 <- bmr$score(predictions = TRUE, ids = TRUE)[x*3-1]$prediction_test[[1]] |>
    data.table::as.data.table()
  df3 <- bmr$score(predictions = TRUE, ids = TRUE)[x*3-2]$prediction_test[[1]] |>
    data.table::as.data.table()
  
  df<- rbind(df1, df2, df3)
  sim_id<- bmr$score(predictions = TRUE, ids = TRUE)[x*3]$task_id
  
  # get max of x and y
  max_xy <- round(max(
    max(df$response, na.rm = TRUE),
    max(df$truth, na.rm = TRUE)
  ) / 10) * 10
  
  p<-  df |>
    ggplot() +
    aes(y = response, x = truth) +
    geom_point(col = "#f08a46", alpha = 0.9) +
    #geom_density_2d(aes(col = after_stat(level))) +
    scale_color_viridis_c(direction = -1, option = "mako") +
    guides(alpha = "none", color = "none") +
    geom_abline(slope = 1) +
    coord_fixed(xlim = c(0, 8), ylim = c(0, 8)) +
    labs(
      title = paste0(sim_id),
      x = "Observed Carbon",
      y = "Modelled Carbon"
    ) +
    theme_linedraw()
  return(p)
}
plot_list <- lapply(y, plot_sims)
final_grid_plot <- wrap_plots(plot_list, ncol = 3, nrow = 4)
ggsave("Plots/grid_of_truth_vs_response.png", plot = final_grid_plot, width = 20, height = 18, dpi = 300)

### importance scores for variables in each simulation
feature_importance<- function (x, xy_lim= 8) {
  df1 <- bmr$score(ids = TRUE)$learner[[x*3]]$importance()|>
    data.table::as.data.table( keep.rownames =  TRUE)
  df2 <- bmr$score(ids = TRUE)$learner[[x*3-1]]$importance()|>
    data.table::as.data.table( keep.rownames =  TRUE)
  df3 <- bmr$score(ids = TRUE)$learner[[x*3-2]]$importance()|>
    data.table::as.data.table( keep.rownames =  TRUE)
  
  df<- rbind(df1, df2, df3)
  
  
  summary_df <- df %>%
    group_by(V1) %>%
    summarise(mean_imp = mean(V2, na.rm = TRUE))%>%
    arrange(desc(mean_imp)) %>%
    mutate(V1 = factor(V1, levels = rev(unique(V1))))
  
  
  n_rows <- min(10, nrow(summary_df))  
  
  plot_df<- summary_df[1:n_rows, ]
  sim_id<- bmr$score(predictions = TRUE, ids = TRUE)[x*3]$task_id
  
  
  
  # Plot
  p<- ggplot(plot_df, aes(x = V1, y = mean_imp)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(
      title = paste0(sim_id),
      x = "variables",
      y = "importance scores"
    ) +
    theme_minimal()
  return(p)
}
plot_list <- lapply(y, feature_importance)
importance_grid_plot <- wrap_plots(plot_list, ncol = 3, nrow = 4)
ggsave("Plots/grid_of_importance_scores.png", plot = importance_grid_plot, width = 20, height = 18, dpi = 300)
