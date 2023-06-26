library(moveVis)

data("whitestork_data")

m <- align_move(m, res = 180, unit = "mins")

frames <- frames_spatial(
  m,
  trace_show = TRUE,
  equidistant = FALSE,
  map_service = "osm",
  map_type = "terrain_bg"
)

frames[[200]] #plots a single frame

frames <- frames %>% 
  add_labels(title = "White Storks (Ciconia ciconia) Migration 2018",
             caption = "Trajectory data: Cheng et al. (2019); 
     Fiedler et al. (2013-2019),https://doi.org/10.5441/001/1.ck04mn78
     Map: OpenStreetMap/Stamen; Projection: Geographic, WGS84",
             x = "Longitude", y = "Latitude") %>%
  add_timestamps(type = "label") %>%
  add_progress(colour = "white") %>%
  add_northarrow(colour = "white", position = "bottomleft") %>% 
  add_scalebar(colour = "black", position = "bottomright", distance = 600)

animate_frames(frames, width = 800, height = 800,
               out_file = "x.gif", end_pause = 1)
