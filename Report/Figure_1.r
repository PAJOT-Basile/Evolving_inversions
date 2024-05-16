# Libraries
require(anyLib)
anyLib(c("tidyverse", "ggmap", "sp", "tmap", "ggpubr", "marmap", "ggspatial", "cowplot"))

# Register for the use of the google maps
register_google(key="")

################## European map ########################
# Prepare the sampling site localisation
LAM_square <- c(48.85637, -5.448272,
                48.85637, -4.354781,
                48.145713, -4.354781,
                48.145713, -5.448272) %>% 
  matrix(ncol=2, byrow=TRUE) %>% 
  as.data.frame() %>% 
  dplyr::rename(lat=V1, lon=V2)

LOK_square <- c(59.24443, 10.56395,
                59.24443, 11.65744,
                58.53378, 11.65744,
                58.53378, 10.56395) %>% 
  matrix(ncol=2, byrow=TRUE) %>% 
  as.data.frame() %>% 
  dplyr::rename(lat=V1, lon=V2)
# Load the european map
europe <- getNOAA.bathy(lon1=-7, lon2=15, lat1=44, lat2=60, resolution = 1)
# Plot the map
europe_map <- autoplot.bathy(europe, geom=c("raster", "contour"), color = "black", lty=0) +
  scale_fill_gradientn(limits = c(europe %>% min, 0),
                       colors = c("lightblue4", "lightblue2")) +
  labs(title="", x="", y="") +
  geom_polygon(data=LAM_square, aes(x=lon, y=lat), color="red", fill=NA, lwd=1.2) +
  geom_polygon(data=LOK_square, aes(x=lon, y=lat), color="red", fill=NA, lwd=1.2) +
  theme_void() +
  theme(legend.position = "none",
        text = element_text(size = 30)) +
  annotation_north_arrow(location = "bl", style = north_arrow_nautical, height = unit(2, "cm"),
                         width = unit(2, "cm"))

################## Lampaul-plouarzel map ########################
# Decide where to put the center of the map
LAM_center <- c(lon = -4.771007, lat = 48.463351)
# Load the map
lam_map <- get_map(location = LAM_center, zoom=17, maptype="satellite")

#Prepare the palette to use to show size evolution along transect
size_palette = c("#4e79a7", "grey75", "#f28e2b")
# Make a function to plot the map
Sampling_LAM <- function(pos, breaks = c(-20, 3, 46, 80), size_palette) {
  get.poly <- function(a, b, r1 = 0.00025, r2 = 0.0005*2.5) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) - 4.771304
    xend     <- r2 * cos(th) - 4.771304
    y        <- r1 * sin(th) + 48.46329
    yend     <- r2 * sin(th) + 48.46329
    data.frame(x, y, xend, yend)
  }
  get.poly_left <- function(a, b, r1 = 0.00025, r2 = 0.0005*2.5) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) - 4.772525
    xend     <- r2 * cos(th) - 4.772525
    y        <- r1 * sin(th) + 48.46242
    yend     <- r2 * sin(th) + 48.46242
    data.frame(x, y, xend, yend)
  }
  get.poly_right_after_simple <- function(a, b, r1 = 0.00025, r2 = 0.0005*2.5) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) - 4.77009
    xend     <- r2 * cos(th) - 4.77009
    y        <- r1 * sin(th) + 48.46417
    yend     <- r2 * sin(th) + 48.46417
    data.frame(x, y, xend, yend)
  }
  get.poly_right_end <- function(a, b, r1 = 0.00025, r2 = 0.0005*2.5) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) -4.77055
    xend     <- r2 * cos(th) -4.77055
    y        <- r1 * sin(th) + 48.46275
    yend     <- r2 * sin(th) + 48.46275
    data.frame(x, y, xend, yend)
  }
  
  df <- rbind(get.poly_left(80, 120) %>% t %>% as.data.frame %>% rev %>% t,
              get.poly(-20, 80),
              get.poly_right_after_simple(-40, -20) %>% arrange(x),
              get.poly_right_end(60, 100)) %>% 
    mutate(coloring=seq_along(1:nrow(.)))
  ggmap(lam_map, base_layer = ggplot(data=df, aes(x=x, y=y))) + 
    geom_segment(aes(xend = xend, yend = yend, color = coloring)) +
    scale_color_gradientn(colors = size_palette,
                          labels = c("8.5", "10.2", "13.1"),
                          breaks = c(nrow(df)/3, 2 * nrow(df)/3, nrow(df)))+
    labs(title = "",
         y="",
         x="",
         tag = "(B)") +
    annotate("text",
             x = -4.7685,
             y = 48.464,
             label = "Abr",
             color = "white",
             size = 10) +
    annotate("text",
             x = -4.7735,
             y = 48.463,
             label = "Exp",
             color = "white",
             size = 10) +
    theme_void() +
    theme(legend.position = "none",
          text = element_text(size = 30))
}
# PLot the map
LAM_map <- Sampling_LAM(pos = 10, size_palette=size_palette)



################## Lökholmen map ########################
# Decide where to put the center of the map
LOK_centre <- c(lat=58.889106, lon=11.110695)
# Load the map
LOK_map <- get_map(LOK_centre, maptype = "satellite", zoom = 17)
# Plot the map
ggmap(LOK_map) +
  labs(title = "Map of sampling location: Lökholmen",
       y="Longitude",
       x="Lattitude") +
  theme(text = element_text(size=20))
# Make a function to load the map
Sampling_LOK <- function(map, breaks = c(-20, 3, 46, 80), size_palette) {
  get.poly_sw <- function(a, b, r1 = 0.0035, r2 = 0.00375) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.1125
    xend     <- r2 * cos(th) + 11.1125
    y        <- r1 * sin(th) + 58.89005
    yend     <- r2 * sin(th) + 58.89005
    data.frame(x, y, xend, yend)
  }
  get.poly_sw_curve <- function(a, b, r1 = 0.000035, r2 = 0.000285) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.10912
    xend     <- r2 * cos(th) + 11.10912
    y        <- r1 * sin(th) + 58.8893
    yend     <- r2 * sin(th) + 58.8893
    data.frame(x, y, xend, yend)
  }
  get.poly_sw_end_curve <- function(a, b, r1 = 0.00075, r2 = 0.001) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.10935
    xend     <- r2 * cos(th) + 11.10935
    y        <- r1 * sin(th) + 58.888623
    yend     <- r2 * sin(th) + 58.888623
    data.frame(x, y, xend, yend)
  }
  get.poly_sw_center_curve <- function(a, b, r1 = 0.000035, r2 = 0.000285) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.10947
    xend     <- r2 * cos(th) + 11.10947
    y        <- r1 * sin(th) + 58.88933
    yend     <- r2 * sin(th) + 58.88933
    data.frame(x, y, xend, yend)
  }
  get.poly_center_curve <- function(a, b, r1 = 0.00035, r2 = 0.0006) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.10938
    xend     <- r2 * cos(th) + 11.10938
    y        <- r1 * sin(th) + 58.889025
    yend     <- r2 * sin(th) + 58.889025
    data.frame(x, y, xend, yend)
  }
  get.poly_center_under_curve <- function(a, b, r1 = 0.00075, r2 = 0.001) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.110205
    xend     <- r2 * cos(th) + 11.110205
    y        <- r1 * sin(th) + 58.8901
    yend     <- r2 * sin(th) + 58.8901
    data.frame(x, y, xend, yend)
  }
  get.poly_center_flat <- function(a, b, r1 = 0.0035, r2 = 0.00375) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.10929
    xend     <- r2 * cos(th) + 11.10929
    y        <- r1 * sin(th) + 58.892698
    yend     <- r2 * sin(th) + 58.892698
    data.frame(x, y, xend, yend)
  }
  get.poly_end <- function(a, b, r1 = 0.000035, r2 = 0.000285) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.11077
    xend     <- r2 * cos(th) + 11.11077
    y        <- r1 * sin(th) + 58.88957
    yend     <- r2 * sin(th) + 58.88957
    data.frame(x, y, xend, yend)
  }
  df <- rbind(
    
    get.poly_sw(-7, -10) %>% arrange(-x),
    get.poly_sw_curve(-7, 40),
    get.poly_sw_end_curve(40, 55),
    get.poly_sw_center_curve(55, 60),
    get.poly_center_curve(60, 72),
    get.poly_center_under_curve(-29.5, -60),
    get.poly_center_flat(-60, -64),
    get.poly_end(-62, -120)
  ) %>% 
    mutate(coloring=seq_along(1:nrow(.)))
  ggmap(map, base_layer = ggplot(data=df, aes(x=x, y=y))) + 
    geom_segment(aes(xend = xend, yend = yend, color = -coloring)) +
    scale_color_gradientn(colors = size_palette,
                          labels = c(8.5, 10.2, 13.1),
                          breaks = c(nrow(df)/3, 2 * nrow(df)/3, nrow(df))) +
    labs(x = "",
         y = "",
         tag = "(A)") + 
    annotate("text",
             x = 11.108,
             y = 58.88925,
             label = "Exp",
             color = "white",
             size = 10) +
    annotate("text",
             x = 11.1125,
             y = 58.88975,
             label = "Abr",
             color = "white",
             size = 10) +
    theme_void() +
    theme(legend.position = "none",
          text = element_text(size = 30))
}
# PLot the map
Map_LOK <- Sampling_LOK(LOK_map, size_palette=size_palette)

################## Merge all maps ########################

ggdraw() +
  coord_equal(xlim = c(-30, 30), ylim = c(-20, 20), expand = FALSE) +
  # Put the europe map
  annotation_custom(ggplotGrob(europe_map), xmin = -5, xmax = 25, ymin = -20, ymax = 20) +
  # Add the Lökholmen map
  annotation_custom(ggplotGrob(Map_LOK), xmin = -30, xmax = -5, ymin = 0, ymax = 15) +
  geom_segment(aes(x = -9.75, xend = 19.35, y = 13.6, yend = 13.25), color = "red", linewidth = 1.3) +
  geom_segment(aes(x = -9.75, xend = 19.35, y = 0.07, yend = 11.8), color = "red", linewidth = 1.3) +
  geom_polygon(data = c(-23.2, 0, -23.2, 13.6, -9.74, 13.6, -9.74, 0) %>% 
                 matrix(ncol = 2, byrow=TRUE) %>% 
                 as.data.frame %>% 
                 rename(x = V1, y = V2), aes(x = x, y = y), color = "red", fill = NA, linewidth = 1.3) +
  # Add the Lampaul map
  annotation_custom(ggplotGrob(LAM_map), xmin = -32.5, xmax = -2.5, ymin = -18.5, ymax = -1.5) +
  geom_segment(aes(x = -9.5, xend = -0.75, y = -4.6, yend = -8.05), color = "red", linewidth = 1.3) +
  geom_segment(aes(x = -9.5, xend = -0.75, y = -18.4, yend = -9.4), color = "red", linewidth = 1.3) +
  geom_polygon(data = c(-9.5, -4.6, -9.5, -18.5, -23.3, -18.5, -23.3, -4.6) %>% 
                 matrix(ncol = 2, byrow=TRUE) %>% 
                 as.data.frame %>% 
                 rename(x = V1, y = V2), aes(x = x, y = y), color = "red", fill = NA, linewidth = 1.3)

