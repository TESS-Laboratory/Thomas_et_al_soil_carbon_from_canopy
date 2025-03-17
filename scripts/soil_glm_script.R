#### script to pick out points for analysis and model them 

#### set up environment ####
#install.packages("performance")
#install.packages("GGally")
#install.packages("fastDummies")
#install.packages("vegan")
#install.packages("mlr3")
#install.packages("mlr3viz")
#install.packages("mlr3learners")
#install.packages("ranger")
#install.packages("mlr3filters")
#install.packages("FSelectorRcpp")
#install.packages("mlr3fselect")
#install.packages("mlr3pipelines")
#install.packages("praznik")
#install.packages("mlr3tuningspaces")
library(scales)
library(terra)
library(ggplot2)
library(GLMMRR)
library(dplyr)
library(performance)
library(GGally)
library(sf)
library(tidyverse)
library(fastDummies)
library(vegan)
library(mlr3)
library(mlr3viz)
library(mlr3learners)
library(ranger)
library(mlr3filters)
library(FSelectorRcpp)
library(mlr3fselect)
library(mlr3filters)
library(mlr3pipelines)
library(praznik)
library(mlr3tuningspaces)


set.seed(42)  # For reproducibility

samples_metrics<- read_sf("Data/soil_samples_w_complete_metrics.fgb")

retained_vars<- samples_metrics%>% select(wmean_cec_5.15cm_mean_1, wmean_cover_z0_5m_2, wmean_classification_2023_1, wmean_percC_5)
####prelim analysis #####
# Create a formula for the GLM where the response column is modeled by all other columns
pairplot <- GGally::ggpairs(retained_vars,  cardinality_threshold = 50)
primary_var <- "wmean_percC_5"
pvar_pos <- match(primary_var, pairplot$xAxisLabels)
plots <- lapply(1:4, function(j){ getPlot(pairplot, i = pvar_pos, j = 1)})
ggmatrix(
  plots,
  nrow = 10,
  ncol = 10,
  xAxisLabels = pairplot$xAxisLabels,
  yAxisLabels = primary_var
)
print(pairplot)
ggsave("var132_140.png", pairplot, path = "Data/Outputs/Pairwaise_plots", width = 35, height = 30, units = "cm")
getPlot(pairplot, i= 3, j= 1)
p_ <- GGally::print_if_interactive


# Create a formula for the GLM where the response column is modeled by all other columns
# Filter out columns with only one unique level (constant columns)


#### PCA Analysis ####



## plot without categorical data 
rda.out <- vegan::rda(retained_vars[,-c(5)], scale = TRUE)
# add scores()
rda_scores <- scores(rda.out)
# add biplot()
biplot(rda.out, type = "text")
## group by depth ### this doesnt work yet
ordihull(rda.out,
         group = Fmetrics_num$Profundidade_cm_,
         col = 1:11,
         lty = 1:11,
         lwd = c(3,6), 
         label = TRUE)


