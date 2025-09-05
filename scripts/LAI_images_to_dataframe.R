####################### script to look at hemishphere R package 
### make this script better for calling locations and adding correct sub names
####################### Environemnt set up #########################
#install.packages("hemispheR")
library(hemispheR)
#install.packages("dplyr")
library(dplyr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages('vegan')
library(vegan)
library(mgcv)
#install.packages("sf")
library(sf)
#install.packages('gstat')
library(gstat)
library(sp)
#library(raster)

################### move data to correct folder ################################
fp<- "C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data"
source_folder <- file.path(fp, "best_LAI")  # Replace with the actual file name

# Define the destination folder path
destination_path <- "C:/Users/jpt215/AppData/Local/R/win-library/4.2/hemispheR/extdata"  # Replace with the actual folder path
# Use file.copy() to copy the folder and its contents
if (dir.exists(source_folder)) {
  # Create the destination folder if it doesn't exist
  if (!dir.exists(destination_path)) {
    dir.create(destination_path, recursive = TRUE)
  }
  
  # Copy the folder and its contents
  if (file.copy(source_folder, destination_path, recursive = TRUE)) {
    cat("Folder copied successfully.\n")
  } else {
    cat("Folder copy failed.\n")
  }
} else {
  cat("Source folder does not exist.\n")
}

############## make uniform light dataframe #############################
#find the camera you need
#list.cameras
# Set the directory to the 'LAI_photo' folder
folder_path <- file.path(source_folder, "Cloudy_LAI")

# List all files in the folder
file_list <- list.files(path = folder_path, full.names = TRUE)
#set up empty dataframe
num_rows<- length(file_list)
Image.data<- data.frame(matrix(nrow = num_rows))


# Loop through each file and read it, analyse, add each value to dataframe
x<- 0
for (file_path in file_list) {
  
  c.im<-system.file(paste0('extdata/', file_path) ,package='hemispheR')
  x<- x+1
  mask.img<- import_fisheye(c.im, circ.mask=camera_fisheye('EOS30D+Sigma-4.5'))
  bin.img<- binarize_fisheye(mask.img)
  gap.img<- gapfrac_fisheye(bin.img, lens='Sigma-4.5',display=TRUE)
  can.tib<- canopy_fisheye(gap.img)
  
  
  # Store the data in a list
  Image.data$ID[x]<- can.tib$id[1]
  Image.data$effective.LAI[x]<- can.tib$Le[1]
  Image.data$actual.LAI[x]<-can.tib$L[1]
  Image.data$clumping.index[x]<- can.tib$LX[1]
  Image.data$LXG1[x]<- can.tib$LXG1[1]
  Image.data$LXG2[x]<- can.tib$LXG2[1]
  Image.data$MTA[x]<- can.tib$MTA.ell[1]
  Image.data$canopy.openness[x] <- can.tib$DIFN[1]
  
}
# clean dataframe
Image.data.clean<- Image.data %>%
  select(-matrix.nrow...num_rows.)%>%
  mutate(light.conditions = "uniform and diffuse")%>%
  mutate(ID = substr(ID,1, nchar(ID) - 4))


##################### make patchy light dataframe #######################
folder_path <- file.path(source_folder, "Sunny_LAI")
# List all files in the folder
file_list <- list.files(path = folder_path, full.names = TRUE)
#set up empty dataframe
num_rows<- length(file_list)
CImage.data<- data.frame(matrix(nrow = num_rows))


# Loop through each file and read it, analyse, add each value to dataframe
x<- 0
for (file_path in file_list) {
  
  c.im<-system.file(paste0('extdata/', file_path) ,package='hemispheR')
  x<- x+1
  mask.img<- import_fisheye(c.im, circ.mask=camera_fisheye('EOS30D+Sigma-4.5'))
  bin.img<- binarize_fisheye(mask.img)
  gap.img<- gapfrac_fisheye(bin.img, lens='Sigma-4.5',display=TRUE)
  can.tib<- canopy_fisheye(gap.img)
  
  
  # Store the data in a list
  CImage.data$ID[x]<- can.tib$id[1]
  CImage.data$effective.LAI[x]<- can.tib$Le[1]
  CImage.data$actual.LAI[x]<-can.tib$L[1]
  CImage.data$clumping.index[x]<- can.tib$LX[1]
  CImage.data$LXG1[x]<- can.tib$LXG1[1]
  CImage.data$LXG2[x]<- can.tib$LXG2[1]
  CImage.data$MTA[x]<- can.tib$MTA.ell[1]
  CImage.data$canopy.openness[x] <- can.tib$DIFN[1]
  
}

# clean dataframe
CImage.data.clean<- CImage.data %>%
  select(-matrix.nrow...num_rows.)%>%
  mutate(light.conditions = "patchy")%>%
  mutate(ID = substr(ID,1, nchar(ID) - 4))


#################### Merge dataframes and export a copy #######################
combined_table <- rbind(Image.data.clean, CImage.data.clean)
colnames(combined_table) <- paste0(colnames(combined_table), "_4")
colnames(combined_table$GEDI.footprints_4) <- paste0("GEDI.footprints")
combined_table<- combined_table %>%
  mutate(ID = gsub("X", "",ID))%>%
  mutate(GEDI.footprints = gsub("GEDI", "",ID))%>%
  mutate(GEDI.footprints = sub("\\..*", "", GEDI.footprints))

#export
write.csv(combined_table, file = file.path(fp, "LAI_table.csv"))

################### combine with plot metadata #####################
#import meta data
combined_table<- read.csv(file.path(fp, "LAI_table.csv"))
GEDI_meta_data<- read.csv(file.path(fp, 'Field_GEDI_plots.csv'))
GEDI_meta_data$GEDI.footprints[!is.na(GEDI_meta_data$Plot)] <- GEDI_meta_data$Plot[!is.na(GEDI_meta_data$Plot)]
#remove sensors because they aren't included here
#make both intergers so you can match the columns up
GEDI_LAI <- combined_table %>%
  filter(GEDI.footprints != "sensor")%>%
  mutate(GEDI.footprints = as.integer(GEDI.footprints))
GEDI_LAI$GEDI.footprints[GEDI_LAI$GEDI.footprints == 100] <- 25
GEDI_LAI$GEDI.footprints[GEDI_LAI$GEDI.footprints == 101] <- 29

#combine
TGEDI_LAI <- GEDI_LAI %>%
  left_join(GEDI_meta_data, by = "GEDI.footprints", multiple = "all")
#remove unwanted columns
TGEDI_LAI<- TGEDI_LAI%>%
  select(-number)%>%
  select(-Date)%>%
  select(-Distinguishing.features..local.info)%>%
  select(-Comments)%>%
  select(-X)%>%
  select(-X.1)%>%
  select(-Lat.Lon..GEDI.)%>%
  select(-Plot)



##read and clean plot location datasets
LAI_plot_locations<- read.csv(file.path (fp, "LAI_locations.csv"))
#LAI_plot_locations$Plot.ID[LAI_plot_locations$Plot.ID >=234 & LAI_plot_locations$Plot.ID<=236] <- 29
height_plot_locations<- read.csv(file.path( fp, "GEDI_height_footprints.csv"))
#height_plot_locations$Plot.ID[height_plot_locations$Plot.ID >=234 & height_plot_locations$Plot.ID<=236] <- 29
height_plot_locations$Plot.ID[height_plot_locations$Plot.ID==233] <- 25
height_plot_locations<- height_plot_locations%>%
  select(Plot.ID, x, y)
LAI_plot_locations<- LAI_plot_locations%>%
  select(Plot.ID, x, y)
plot_locations<- rbind(LAI_plot_locations, height_plot_locations)
## fix field mistake and add plot 29
plot_29<- c(29, -64.390240000000006, -12.104649999999999)
plot_locations<- rbind(plot_locations, plot_29)
#create a list of all the plot names
plot_locations$Plot.ID<- as.integer(plot_locations$Plot.ID)
plot_list<- unique(plot_locations$Plot.ID)
#create a function that takes every value in a list and analyses it using a dataframe
subsetAndSaveTables <- function(dataframe, values_list) {
  for (value in values_list) {
    
    # subset a column in the dataframe called 'Plot' for any rows that match the current value in the list
    subset_df <- dataframe[dataframe$Plot.ID == value, ]    
    #aveage the subsetted locations
    xmean <- mean(subset_df$x)
    
    ymean<- mean(subset_df$y)
    # add average values to dataframe
    TGEDI_LAI$plot.x[TGEDI_LAI$GEDI.footprints== value]<- xmean
    TGEDI_LAI$plot.y[TGEDI_LAI$GEDI.footprints== value]<- ymean
    assign('coord_table', TGEDI_LAI, envir = .GlobalEnv)
  }
}

subsetAndSaveTables(plot_locations, plot_list)
#fix degredation state column
coord_table<- coord_table%>%
  mutate(Degradation = ifelse(Degradation == 'unkown', 'Structurally intact', Degradation))%>%
  mutate(Degradation = ifelse(Degradation == 'burned', 'Burned', Degradation))%>%
  mutate(Degradation = ifelse(Degradation == 'logged', 'Logged', Degradation))
coord_table
#export
write.csv(coord_table, file = "LAI_w_meta_data.csv")