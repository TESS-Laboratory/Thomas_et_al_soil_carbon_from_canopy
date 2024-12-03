####################### script to look at hemishphere R package 

####################### Environemnt set up #########################
install.packages("hemispheR")
library(hemispheR)
install.packages("dplyr")
library(dplyr)
install.packages("ggplot2")
library(ggplot2)
install.packages('vegan')
library(vegan)
library(mgcv)
install.packages("sf")
library(sf)
install.packages('gstat')
library(gstat)
library(sp)
library(raster)

################### move data to correct folder ################################
# Define the source file name (in the working directory)
source_folder <- "best_LAI"  # Replace with the actual file name

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
folder_path <- "best_LAI/Cloudy_LAI"

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
folder_path <- "best_LAI/Sunny_LAI"

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
combined_table<- combined_table %>%
  mutate(ID = gsub("X", "",ID))%>%
  mutate(GEDI.footprints = gsub("GEDI", "",ID))%>%
  mutate(GEDI.footprints = sub("\\..*", "", GEDI.footprints))

#export
write.csv(combined_table, file = "LAI_table.csv")

################### combine with plot metadata #####################
#import meta data
combined_table<- read.csv("LAI_table.csv")
GEDI_meta_data<- read.csv('Field_GEDI_plots.csv')
GEDI_meta_data$GEDI.footprints[!is.na(GEDI_meta_data$Plot)] <- GEDI_meta_data$Plot[!is.na(GEDI_meta_data$Plot)]
#remove sensors because they aren't included here
#make both intergers so you can match the columns up
GEDI_LAI <- combined_table %>%
  filter(GEDI.footprints != "sensor")%>%
  mutate(GEDI.footprints = as.integer(GEDI.footprints))
GEDI_LAI$GEDI.footprints[GEDI_LAI$GEDI.footprints == 100] <- 25
GEDI_LAI$GEDI.footprints[GEDI_LAI$GEDI.footprints == 101] <- 29

#comine
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
LAI_plot_locations<- read.csv("LAI_locations.csv")
#LAI_plot_locations$Plot.ID[LAI_plot_locations$Plot.ID >=234 & LAI_plot_locations$Plot.ID<=236] <- 29
height_plot_locations<- read.csv("GEDI_height_footprints.csv")
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

################# finding distance from edge ################
##load in shapefile
RCsite<- st_read("RIO_CAUTARIO.shp")
RCbound<- st_boundary(RCsite)

comp.coord<- coord_table[complete.cases(coord_table[, c("plot.x", "plot.y")]), ]
xy_sf <- st_as_sf(comp.coord, coords = c("plot.x", "plot.y"), crs = st_crs(RCbound))
# Calculate minimum distances
distances <- st_distance(xy_sf, RCbound)

# Add distances to the original table
comp.coord$min_distance <- distances

Ccomp.coord<- comp.coord%>%
  select(-X.2,-X.1, -X, -Degradation.....unknown., -Degradation.....unkown., -X.Structurally.Intact., -X.Structurally.Intact..1)
# Optionally, save the modified table
write.csv(Ccomp.coord, "LAI_w_meta_data.csv")


############# LAI interpolation and comparison ##################

# create a data frame called 'known_data' with columns 'LAI', 'x', and 'y'
known_data<- coord_table%>%
  filter(State == "Intact")%>%
  select(actual.LAI, plot.x, plot.y)
# 'LAI' represents the LAI values, and 'x' and 'y' are the coordinates

# Create a spatial points data frame
coordinates(known_data) <- ~plot.x + plot.y

# Define the grid of points where you want to predict LAI values
secPoints<- coord_table%>%
  filter((State == "Secondary"))%>%
  dplyr::select(plot.x, plot.y, ID, actual.LAI, Age)

# Replace 'grid_x' and 'grid_y' with the desired grid points

grid_x <- seq(min(secPoints$plot.x), max(secPoints$plot.x), length.out = 100)
grid_y <- seq(min(secPoints$plot.y), max(secPoints$plot.y), length.out = 100)
grid_points <- expand.grid(x = grid_x, y = grid_y)
coordinates(grid_points) <- ~x + y
# Perform IDW interpolation
idw_model <- idw(actual.LAI ~ 1, known_data, grid_points)

# Predict LAI values for the grid points and add them to a column
grid_points$LAI_pred <- idw_model$var1.pred

#convert back to a dataframe
grid_points_df <- as.data.frame(grid_points)

# You can plot the interpolated values
ggplot(grid_points_df, aes(x, y, fill = LAI_pred)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(fill = "Predicted LAI") +
  theme_minimal()

#create a raster to extract values from
extent <- extent(grid_points)

# Define the number of rows and columns for the raster grid
nrows <- 100  # Adjust as needed
ncols <- 100  # Adjust as needed

# Create an empty raster with the specified extent and grid size
raster_grid <- raster(extent, nrow = nrows, ncol = ncols)
# populate
populated_raster <- rasterize(grid_points, raster_grid, field = "LAI_pred")
plot(populated_raster)

new_coords <- data.frame(x = secPoints$plot.x, y = secPoints$plot.y)

# Create a spatial points data frame for the new coordinates
coordinates(new_coords) <- ~x + y

# Extract predicted LAI values and add to sec points table
secPoints$predicted_LAI <- extract(populated_raster, new_coords)

## calculate difference between predicted and observed LAI
secPoints$LAI.difference<- secPoints$actual.LAI - secPoints$predicted_LAI

##plot how that difference changes with age 
plot(secPoints$Age, secPoints$LAI.difference)
LAI.diffGAM<- gam(LAI.difference ~ s(Age, k= ), data = secPoints, method = "REML")
summary(LAI.diffGAM)
gam.check(LAI.diffGAM)
concurvity(LAI.diffGAM, full = FALSE)
plot(LAI.diffGAM, pages = 1, residuals = TRUE, pch = 1, cex = 1, shade= TRUE)

################### Quick comparisons ########################
coord_table<- read.csv('LAI_w_meta_data.csv')

## effect of light conditions
ggplot(data = coord_table, aes(x = light.conditions, y = actual.LAI, fill = light.conditions)) +
  geom_boxplot() +
  labs(x = "Degradation", y = "LAI Value", title = "Box Plot of LAI by Degradation") +
  theme_minimal()
light.lai<- coord_table$actual.LAI[coord_table$light.conditions == 'patchy']
cloud.lai<- coord_table$actual.LAI[coord_table$light.conditions == 'uniform and diffuse']
t.test(light.lai, cloud.lai)

##LAI distribution
hist(coord_table$actual.LAI)
shapiro.test(coord_table$actual.LAI)
ks.test(coord_table$actual.LAI, "pnorm", mean = mean(coord_table$actual.LAI), sd = sd(coord_table$actual.LAI))

#############LAI vs age ##################
coord_table$Age<- as.numeric(coord_table$Age)
coord_table$actual.LAI<- as.numeric(coord_table$actual.LAI)

#########this tested for linear regression which does not work 
#plot(coord_table$Age,coord_table$actual.LAI)
#age.lAi.m<- lm(actual.LAI~ Age, data = coord_table, na.action = na.exclude)
#abline(age.lAi.m)
#summary(age.lAi.m)
#cor.test(coord_table$Age, coord_table$actual.LAI, method = 'kendall')

### polynomial regression
coord_table_clean <- coord_table[complete.cases(coord_table[, c("Age", "actual.LAI")]), ]
fit <- lm(actual.LAI ~ poly(Age, 2), data = na.omit(coord_table_clean))
summary(fit)

#### plot with confidence intervals
# Make predictions using the model
coord_table_clean$Predicted_LAI <- predict(fit)

# Create the plot
ggplot(coord_table_clean, aes(x = Age, y = actual.LAI)) +
  geom_point() +  # Plot the data points
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "blue", se = TRUE) +  # Plot the model with confidence intervals
  labs(title = "Polynomial Regression with Confidence Intervals") +
  theme_minimal()

##### look at how this model works for a 50 year old tree
x_new <- data.frame(Age = 50)
predicted_y <- predict(fit, newdata = x_new)

# Print the predicted value
cat("Predicted y for x = 50:", predicted_y)
##### it aint right. thats for sure. 

############### LAI vs degredation   #####################
### boxplots of LAI for various degredation types
ggplot(data = coord_table, aes(x = Degradation, y = actual.LAI, fill = Degradation)) +
  geom_boxplot() +
  labs(x = "Degradation", y = "LAI Value", title = "Box Plot of LAI by Degradation") +
  theme_minimal()

###T_tests to see if the LAI of degredation types is significantly different
burned.lai<- coord_table$actual.LAI[coord_table$Degradation == 'Burned']
logged.lai<- coord_table$actual.LAI[coord_table$Degradation == 'Logged']
intact.lai<- coord_table$actual.LAI[coord_table$Degradation == 'Structurally intact']
t.test(intact.lai, burned.lai, alternative = 'greater')
t.test(intact.lai, logged.lai, alternative = 'greater')
t.test(burned.lai, logged.lai, alternative = 'less')

###########LAI vs distance from edge ################
plot(coord_table$min_distance, coord_table$actual.LAI)
distlm<- lm(coord_table$actual.LAI ~ coord_table$min_distance)
abline(distlm)
summary(distlm)
########## GAM for LAI #################
## change character columns to factors
coord_table <- coord_table%>%
  #select(-X.1, -X.2 ) %>%               # Remove the "X1" column
  mutate(across(c("light.conditions", "Degradation", "State"), as.factor)) %>%
  mutate(Age_rectified = as.numeric(Age_rectified))%>%
  mutate(Age= as.numeric(Age))
#coord_table_clean <- coord_table[complete.cases(coord_table[, c("Age_rectified", "actual.LAI", "Degradation", "State")]), ]

### set up GAM with all independent variables
big.gam<- gam(actual.LAI ~ s(min_distance) + s(Age) + Degradation, data = coord_table, method = "REML")
summary(big.gam)
gam.check(big.gam)
concurvity(big.gam, full = FALSE)
plot(big.gam, pages = 1, residuals = TRUE, pch = 1, cex = 1, shade= TRUE)



###effect size
# Create a new dataset for plotting the effect
new_data <- data.frame(Age = seq(min(coord_table$Age), max(coord_table$Age)),
                       Degradation = "Logged")  # Specify a category
                       

# Predict the effect using the GAM model
predictions <- predict(Sp.gam, new_data, type = "link", se.fit = TRUE)

# Create a plot of the effect with confidence intervals
plot(new_data$Age, exp(predictions$fit), type = "l", lwd = 2, 
     xlab = "Age", ylab = "Effect on actual.LAI", 
     ylim = c(0, max(exp(predictions$fit + 1.96 * predictions$se.fit))),
     main = "Effect of Age on actual.LAI")

# Add confidence intervals
lines(new_data$Age, exp(predictions$fit + 1.96 * predictions$se.fit), lty = 2)
lines(new_data$Age, exp(predictions$fit - 1.96 * predictions$se.fit), lty = 2)
# make table for DCA
DCA_table<- GEDI_LAI%>%
  select(actual.LAI, GEDI.footprints, MTA, canopy.openness, LXG1)%>%
  na.omit()

DCAanalysis<-decorana(DCA_table)
#plot
plot(DCAanalysis,display="species")
##shows that LXG1 and MTA arent having much of an effect but canopy.openess and actual.LAI work strongly in the same vector but opposite directions


######################## error testing ####################
c.im<-system.file('extdata/best_LAI/Cloudy_LAI/2.3.jpg',package='hemispheR')
c.im |>
  import_fisheye(circ.mask=camera_fisheye('EOS30D+Sigma-4.5'))|>
  binarize_fisheye()|>
  gapfrac_fisheye(lens='Sigma-4.5',display=TRUE)|>
  canopy_fisheye()

c.im<-system.file(file_path ,package='hemispheR')
mask.img<- import_fisheye(c.im, circ.mask=camera_fisheye('EOS30D+Sigma-4.5'))
bin.img<- binarize_fisheye(mask.img)
gap.img<- gapfrac_fisheye(bin.img, lens='Sigma-4.5',display=TRUE)
can.tib<- canopy_fisheye(gap.img)



TGEDI_LAI$plot_x[TGEDI_LAI$GEDI.footprints== value]<- xmean
TGEDI_LAI$plot.y[TGEDI_LAI$GEDI.footprints== value]<- ymean

##not a useful plot
ggplot(data = GEDI_LAI, aes(x = 1:length(actual.LAI), y = actual.LAI, color = Degradation)) +
  geom_point() +
  scale_color_manual(values = c("logged" = "blue", "burned" = "green", "NA" = "red")) +
  labs(x = "Index", y = "LAI Value", title = "LAI vs. Degradation") +
  theme_minimal()
### Spatial GAM
Sp.gam<- gam(actual.LAI ~ s(plot_x, plot.y, k=4)+ s(Age)+ Degradation, data = coord_table, method = 'REML')
summary(Sp.gam)
gam.check(Sp.gam)
plot(Sp.gam, scheme = 2)