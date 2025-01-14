set.seed(42)

##drop useless columns
non_sf<-samples_metrics%>%
  st_drop_geometry()%>%
  select(c(-"Codigo", -"n", -"GEDI.footprints",-"id_clean" ))

## make everything numeric or factor
convert_to_numeric_or_factor <- function(df) {
  df %>%
    mutate(across(everything(), ~ {
      if (suppressWarnings(all(!is.na(as.numeric(.))))) {
        as.numeric(.)
      } else {
        as.factor(.)
      }
    }))
}
converted_data <- convert_to_numeric_or_factor(non_sf)
##make sure headers work in mlr3 ecosystem
colnames(converted_data)<-make.names(colnames(converted_data))

###sort out incomplete rows

##Create a task for each data type level  
tsk_all_data = as_task_regr(converted_data, target = "percC",
                               id = "soils")
tsk_all_data$set_col_roles("Profundidade_cm_", roles = "group")
RS_data_subset<- converted_data%>% select(-"x15N",-"percN",-"x13C",-"CperN", -"effective.LAI",-"actual.LAI", -"clumping.index",
                                           -"LXG1",-"LXG2",-"MTA",-"canopy.openness",-"light.conditions")
tsk_RS_data = as_task_regr(RS_data_subset, target = "percC",
                            id = "soils")
tsk_RS_data$set_col_roles("Profundidade_cm_", roles = "group")
FA_data_subset<- converted_data%>% select("bdod_0.5cm_mean", "bdod_5.15cm_mean","bdod_15.30cm_mean",      
                                          "bdod_30.60cm_mean","bdod_60.100cm_mean","bdod_100.200cm_mean","cec_0.5cm_mean",         
                                          "cec_5.15cm_mean","cec_15.30cm_mean","cec_30.60cm_mean","cec_60.100cm_mean",      
                                          "cec_100.200cm_mean","clay_0.5cm_mean","clay_5.15cm_mean","clay_15.30cm_mean",     
                                          "clay_30.60cm_mean","clay_60.100cm_mean","clay_100.200cm_mean","nitrogen_0.5cm_mean",    
                                          "nitrogen_5.15cm_mean", "nitrogen_15.30cm_mean","nitrogen_30.60cm_mean","nitrogen_60.100cm_mean", 
                                          "nitrogen_100.200cm_mean", "ocs_0.30cm_mean","phh2o_0.5cm_mean","phh2o_5.15cm_mean",     
                                          "phh2o_15.30cm_mean","phh2o_30.60cm_mean","phh2o_60.100cm_mean","phh2o_100.200cm_mean",  
                                          "sand_0.5cm_mean","sand_5.15cm_mean" ,       "sand_15.30cm_mean"   ,    "sand_30.60cm_mean",      
                                          "sand_60.100cm_mean"   ,   "sand_100.200cm_mean"   ,  "soc_0.5cm_mean"  ,        "soc_5.15cm_mean" ,       
                                          "soc_15.30cm_mean"  ,      "soc_30.60cm_mean"    ,    "soc_60.100cm_mean"   ,    "soc_100.200cm_mean",
                                          "HydroRiver_raster", "NDVI_variance",           "Profundidade_cm_","Age" ,"State",
                                          "Degradation" , "plot.x" ,                 "plot.y" ,                 "min_distance", "percC", "Profundidade_cm_")
tsk_FA_data = as_task_regr(FA_data_subset, target = "percC",
                           id = "soils")
tsk_FA_data$set_col_roles("Profundidade_cm_", roles = "group")

## feature selection 
all_flt_gain = flt("information_gain")
RS_flt_gain = flt("information_gain")
FA_flt_gain = flt("information_gain")
all_IG<- all_flt_gain$calculate(tsk_all_data)
as.data.table(all_IG)
RS_IG<-RS_flt_gain$calculate(tsk_RS_data)
as.data.table(RS_IG)
FA_IG<-FA_flt_gain$calculate(tsk_FA_data)
as.data.table(FA_IG)


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
tasks = tsks(c("tsk_all_data", "tsk_RS_data", "tsk_FA_data"))
learners = lrns(c("regr.ranger", "regr.featureless"), predict_type = "response")
rsmp_ho = rsmp("holdout")

design = benchmark_grid(tasks, learners, rsmp_ho)
head(design)
bmr = benchmark(design)
bmr$score()