##### a script for investigating the relationships between microclimates in the forest and forest structure###

###set up environment ####
install.packages("tidyverse")
install.packages("lme4")
install.packages("ggeffects")
install.packages("patchwork")
library(tidyverse)
library(lme4)
library(ggeffects)
library(patchwork)

#### read in data ####
sensor_table<- read_csv("C:/workspace/PhD year 2/datasets/Rio_Cautario/climate stuff/Daily_sensor_spatial_dataset.csv")
# Function to clean column names by removing spaces and brackets
clean_headers <- function(df) {
  clean_names <- gsub("[[:space:]]", "", colnames(df))  # Remove all spaces
  clean_names <- gsub("\\(|\\)", "_", clean_names)  # Remove all brackets
  clean_names<- gsub("%", "perc", clean_names) #remove percent sign
  clean_names<- gsub("/", "per", clean_names)
  clean_names<- gsub("²", "sq", clean_names)
  colnames(df) <- clean_names
  return(df)
}

# Apply the function to the sample data frame
sensor_table <- clean_headers(sensor_table)
##### big model ###
full_model<- glm(temp_degC_max ~ effective.LAI*canopy.openness*DFE..m.*CS_rad*elevation*slope*aspect*Temp.Max._C_*Chuva_mm_*Umi.Min._perc_*Vel.Vento_mpers_*Radiacao_KJpermsq_, data= sensor_table, family = Gamma(link = "log") )

                    