library(moveVis)
library(move)
library(raster)

# loading example data as data.frame 
data("whitestork_data", package = "moveVis")
# includes only what we need for the animation: location, timestamp and individual ID/name
colnames(df)
range(df$timestamp)

# often, individual IDs look like this
unique(df[["individual-local-identifier"]])

# not very readable for a animation legend, so we create a shorter version of these names:
df$name <- sapply(df$`individual-local-identifier`, function(x)
  strsplit(x, " /")[[1]][1], USE.NAMES = F
)
unique(df$name)
df$name <- gsub("-", " ", gsub("[+]", "", gsub(" ", "", df$name)))
unique(df$name)

# now, we are ready to convert the data.frame to a spatial object
# importatn: define the correct CRS in which your coordinates where recorded
m <- df2move(
  df, proj = "+proj=longlat +datum=WGS84", x = "location-long", 
  y = "location-lat", time = "timestamp", track_id = "name", 
  removeDuplicatedTimestamps = TRUE
)

# lets check the sampling rate first:
lag <- unlist(timeLag(m, unit = "mins"))
median(lag)
sd(lag)
# sampling rate is roghly 5 minutes and differing over time with a SD of 36.86 minutes
# indicating coverage gaps = everything is not perfectly regular.
# however: animations use a discrete speed for switching frames.

# moveVis needs to assign each location of a trajectory to a specific frame, 
# the sampling times of all locations across all trajectories need to be aligned 
# to share a uniform temporal resolution and uniform time stamps that can be assigned to frames.

# lets choose a rate that is coarse enough to not get too many frames in the end: every 3 hours
m <- align_move(m, res = 180, unit = "mins")
length(unique(timestamps(m)))

# next, we want to use custom imagery as a dynamically changing basemap
# we use MODIS NDVI, which we load as a temporal stack of rasters
ndvi <- unstack(stack("data/MOD13Q1_NDVI_latlon.tif"))
#  according acquisition times per raster are defined in a separate vector.
ndvi_times <- as.POSIXct(c(
  "2018-07-28", "2018-08-13", "2018-08-29", "2018-09-14"
))
ndvi

# lastly, lets define a manual extent for which we want to animate:
ext <- extent(m)*1.1
ext@xmin <- ext@xmin*1.3
ext@xmax <- ext@xmax*1.3

# crop
#ndvi <- crop(ndvi, ext*1.2)

# now we ar ready to create frames from our movement trajectories and the rasters
frames <- frames_spatial(
  m, r_list = ndvi, r_times = ndvi_times,
  fade_raster = T, ext = ext,
  trace_show = T, trace_colour = "white"
)

frames[[400]]

# map colours can be overwritten to apply a custom colour ramp
frames <- frames %>% 
  add_colourscale("gradient", colours = c("sandybrown", "white",
                                          "green", "darkgreen"), na.colour = "lightskyblue",
                  legend_title = "NDVI")

frames[[200]]

# customize frames:
frames <- frames %>% 
  add_labels(title ="White Storks (Ciconia ciconia) Migration 2018",
             caption = "Trajectory data: Cheng et al. (2019);
    Fiedler et al. (2013-2019), doi:10.5441/001/1.ck04mn78 Map: NASA
    MODIS MOD13Q1 NDVI/Stamen; Projection: Geographic, WGS84",
             x = "Longitude", y = "Latitude") %>%
  add_timestamps(type = "label") %>% 
  add_progress(colour = "white") %>%
  add_northarrow(colour = "black", position = "bottomleft") %>% 
  add_scalebar(colour = "black", position = "bottomright",
               distance = 600)

frames[[200]]

# animate!
animate_frames(frames, width = 800, height = 800,
               out_file = "white_storks_ndvi.mov", end_pause = 1)
