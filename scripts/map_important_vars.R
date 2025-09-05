
### make plots of maps to spatially interrogate understanding 
### import complete raster stack from complire rasters.R
fp<- "C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Soil_manuscript_data"
metrics<- rast(file.path(fp, "combined_metrics_raster.tif"))
study_site <- st_read("C:/Users/jpt215/OneDrive - University of Exeter/PhD_Data/Rio_Cautario/RC_boundary_EPSG4326/RC_boundary_EPSG4326.shp", crs =st_crs("EPSG:4326"))
prj<- crs(metrics)

x_sf <- st_transform(study_site$geometry, crs = st_crs(prj))
#####set up environment
library(terra)
library(rsi)
library(patchwork)
library(ggplot2)
library(viridis)
library(sf)
## make list of variable names you wish to map (here vars from parsimonius GLM were chosen)
variables<- c("isd_3","year_of_last_fire_1","imax_3", "imean_3", "zskew_3","LAD_3","ocs_0-30cm_mean_1", "bdod_15-30cm_mean_1")

map_plot_func<- function(layers, df=metrics){
  layer<- layers[1]
  l<-df[[layer]]
  table<- as.data.frame(l)
  q2<- quantile(table[[1]], probs = 0.99)
  q1<- quantile(table[[1]], probs = 0.01)
  q2<-round(q2)
  q1<-round(q1)
  r_df <- as.data.frame(l, xy = TRUE)  
  names(r_df)[3] <- "value"
  
  p<-ggplot() +
    geom_raster(data = r_df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = x_sf, fill = NA, color = "red", size = 0.5) + # Outline in red
    #ggtitle(layer)+
    scale_fill_viridis(limits = c(q1, q2)) +
    theme_minimal()
    
    
  return(p)
}

plot_list<- lapply(variables, map_plot_func)

map_grid<-wrap_plots(plot_list, ncol = 3, nrow = 3)
map_grid
ggsave("Plots/grid_of_variable_maps.png", plot = map_grid, width = 20, height = 18, dpi = 300)
