
### make plots of maps to spatially interrogate understanding 
### import complete raster stack from complire rasters.R
metrics<- rast("C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Large_Data/combined_metrics_raster.tif")

#####set up environment
library(terra)
library(rsi)
library(patchwork)
library(ggplot2)
library(viridis)
## make list of variable names you wish to map (here vars from parsimonius GLM were chosen)
variables<- c("isd_3","year_of_last_fire_1","imax_3", "imean_3", "zskew_3","LAD_3","ocs_0-30cm_mean_1", "bdod_15-30cm_mean_1")

map_plot_func<- function(layers, df=metrics){
  layer<- layers[1]
  l<-df[[layer]]
  
  r_df <- as.data.frame(l, xy = TRUE)  
  names(r_df)[3] <- "value"
  
  p<-ggplot(r_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_equal() +
    ggtitle(layer)+
    scale_fill_viridis(trans = "log") +
    theme_minimal()
  return(p)
}

plot_list<- lapply(variables, map_plot_func)

map_grid<-wrap_plots(plot_list, ncol = 3, nrow = 3)
map_grid
ggsave("Plots/grid_of_variable_maps.png", plot = map_grid, width = 20, height = 18, dpi = 300)