# Libraries
require(anyLib)
anyLib(c("tidyverse", "ggmap", "sp", "tmap", "ggpubr", "marmap", "ggspatial"))

# Register for the use of the google maps
register_google(key="AIzaSyD8HGsK8KvX2MMO5zIPzgehipwm29Fdvc0")

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
autoplot.bathy(europe, geom="contour", color="grey", lwd=0.1, lty=0) +
  labs(title="Map of sampling sites in Europe: Sweden and Britanny",
       x="", y="") +
  geom_polygon(data=LAM_square, aes(x=lon, y=lat), color="red", fill=NA, lwd=1) +
  geom_polygon(data=LOK_square, aes(x=lon, y=lat), color="red", fill=NA, lwd=1) +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        text = element_text(size=20)) +
  annotation_north_arrow(location = "bl", style = north_arrow_nautical)

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
    x        <- r1 * cos(th) -4.771304
    xend     <- r2 * cos(th) -4.771304
    y        <- r1 * sin(th) + 48.46329
    yend     <- r2 * sin(th) + 48.46329
    data.frame(x, y, xend, yend)
  }
  get.poly_left <- function(a, b, r1 = 0.00025, r2 = 0.0005*2.5) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 500)
    x        <- r1 * cos(th) -4.772525
    xend     <- r2 * cos(th) -4.772525
    y        <- r1 * sin(th) + 48.46242
    yend     <- r2 * sin(th) + 48.46242
    data.frame(x, y, xend, yend)
  }
  get.poly_right_after_simple <- function(a, b, r1 = 0.00025, r2 = 0.0005*2.5) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 250)
    x        <- r1 * cos(th) -4.77009
    xend     <- r2 * cos(th) -4.77009
    y        <- r1 * sin(th) + 48.46417
    yend     <- r2 * sin(th) + 48.46417
    data.frame(x, y, xend, yend)
  }
  get.poly_right_end <- function(a, b, r1 = 0.00025, r2 = 0.0005*2.5) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 333)
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
    geom_segment(aes(xend = xend, yend = yend, color = coloring), alpha=0.5) +
    scale_color_gradientn(colors = size_palette,
                          labels = c("8.5", "10.2", "13.1"),
                          breaks = c(nrow(df)/3, 2 * nrow(df)/3, nrow(df)))+
    labs(title = "Map of sampling location: Lampaul-plouarzel",
         color = "Sell size (mm)",
         y="Longitude",
         x="Lattitude") +
    theme(text=element_text(size=20))
}
# PLot the map
Sampling_LAM(pos = 10, size_palette=size_palette)



################## Lökholmen map ########################
# Decide where to put the center of the map
LOK_centre <- c(lat=58.889106, lon=11.110695)
# Load the map
LOK_map <- get_map(LOK_centre, maptype = "satellite", zoom = 16)
# Plot the map
ggmap(LOK_map) +
  labs(title = "Map of sampling location: Lökholmen",
       y="Longitude",
       x="Lattitude") +
  theme(text = element_text(size=20)) +
  geom_segment(get.poly())
# Make a function to load the map
Sampling_LOK <- function(map, breaks = c(-20, 3, 46, 80), size_palette) {
  get.poly_sw <- function(a, b, r1 = 0.0035, r2 = 0.0037) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.10878
    xend     <- r2 * cos(th) + 11.10878
    y        <- r1 * sin(th) + 58.89145
    yend     <- r2 * sin(th) + 58.89145
    data.frame(x, y, xend, yend)
  }
  get.poly_sw_sc1 <- function(a, b, r1 = 0.00035, r2 = 0.00055) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.110425
    xend     <- r2 * cos(th) + 11.110425
    y        <- r1 * sin(th) + 58.88776
    yend     <- r2 * sin(th) + 58.88776
    data.frame(x, y, xend, yend)
  }
  get.poly_sc1 <- function(a, b, r1 = 0.00035, r2 = 0.00055) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.11065
    xend     <- r2 * cos(th) + 11.11065
    y        <- r1 * sin(th) + 58.88864
    yend     <- r2 * sin(th) + 58.88864
    data.frame(x, y, xend, yend)
  }
  get.poly_sc2 <- function(a, b, r1 = 0.0035, r2 = 0.0037) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.110565
    xend     <- r2 * cos(th) + 11.110565
    y        <- r1 * sin(th) + 58.89179
    yend     <- r2 * sin(th) + 58.89179
    data.frame(x, y, xend, yend)
  }
  get.poly_sc2_sc3 <- function(a, b, r1 = 0.00175, r2 = 0.00195) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.11076
    xend     <- r2 * cos(th) + 11.11076
    y        <- r1 * sin(th) + 58.8900475
    yend     <- r2 * sin(th) + 58.8900475
    data.frame(x, y, xend, yend)
  }
  get.poly_sc3 <- function(a, b, r1 = 0.00035, r2 = 0.00055) {
    th.start <- pi * (1 - a / 100)
    th.end   <- pi * (1 - b / 100)
    th       <- seq(th.start, th.end, length = 1000)
    x        <- r1 * cos(th) + 11.110965
    xend     <- r2 * cos(th) + 11.110965
    y        <- r1 * sin(th) + 58.8886675
    yend     <- r2 * sin(th) + 58.8886675
    data.frame(x, y, xend, yend)
  }
  df <- rbind(
    get.poly_sw(-60, -62.5),
    get.poly_sw_sc1(30, 60),
    get.poly_sc1(-45, -52),
    get.poly_sc2(-51, -53),
    get.poly_sc2_sc3(-53, -55),
    get.poly_sc3(-58, -85)
  ) %>% 
  mutate(coloring=seq_along(1:nrow(.)))
  ggmap(map, base_layer = ggplot(data=df, aes(x=x, y=y))) + 
    geom_segment(aes(xend = xend, yend = yend, color = coloring)) +
    scale_color_gradientn(colors = size_palette,
                          labels = c("8.5", "10.2", "13.1"),
                          breaks = c(nrow(df)/3, 2 * nrow(df)/3, nrow(df)))
}
# PLot the map
Sampling_LOK(LOK_map, size_palette=rev(size_palette))


ggarrange(ggmap(LOK_map) +
  labs(title = "Map of sampling site: Lökholmen")+
  theme(text = element_text(size=20)),
rbind(LOKn, LOKs) %>% 
  ggplot(aes(LCmeanDist, Length, color=Length)) +
  geom_point(size=3, alpha=0.7) +
  scale_colour_gradientn(colors=size_palette, 
                         limits=c(min(data$Length, na.rm=TRUE) - 0.01, 
                                  max(data$Length, na.rm=TRUE) + 0.01),
                         breaks=c(8, 13, 18)) +
  labs(x="Position along the transect (m)",
       y="Shell size (mm)",
       title = "Distribution of shell size along the transect") +
  theme_bw(), nrow=2, heights = c(3, 1), widths=c(2, 1))

