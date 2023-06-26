#
# Example analysis for "Revisitation analysis uncovers spatio-temporal patterns in animal movement data"
# using the recurse package
#
# author: Chloe Bracis

#-------------------------------------
#load libraries
#-------------------------------------
library(recurse)
library(move)
library(ggplot2)
library(ggmap)
library(RgoogleMaps)
library(raster)
library(scales)
library(viridis)
library(lubridate)
library(reshape2)
library(raster)
library(rworldmap)
library(maptools)
library(cluster)

#-------------------------------------
# supporting functions
#-------------------------------------

####
#' Assigns each location/datetime pair to day or night
#'
#' @param locations matrix of (x,y) locations in longiture, latitude
#' @param datetime POSIXct datetime of observation in the local timezone
#'
#' @return vector of day and night 
#'
calculateTimeOfDay = function(locations, datetime)
{
	require(maptools)
	
	sunrise = sunriset(locations, dateTime = datetime, direction = "sunrise", POSIXct.out = TRUE)$time
	sunset = sunriset(locations, dateTime = datetime, direction = "sunset", POSIXct.out = TRUE)$time

	light = ifelse( datetime < sunrise, "night",
					ifelse( datetime < sunset, "day", "night" ) )
	
	return(light)
}

#-------------------------------------
# Download data
#-------------------------------------

# download datafile from Movebank: http://dx.doi.org/10.5441/001/1.46ft1k05
# the data on movebank is 19 vultures, we select just Leo to analyze
vultures = split(move("Turkey vultures in North and South America (data from Dodge et al. 2014).csv"))
leo = vultures$Leo

# the file LeoRoadDistance.csv is available in the online supplementary material
# the National Road Network (NRN) Shapefile for Saskatchwan was obtained from www.geobase.ca
# we projected it to UTM-13, then used the gDistance() function in the rgeos package 
# to calculate the shortest distance from each trajectory location to a road
roadDist = read.csv("LeoRoadDistance.csv")


#-------------------------------------
# Initial plotting, radius selection
#-------------------------------------

# plot leo on a map
leomap = getMap(resolution = "low")
plot(leomap, xlim = range(leo@coords[,1]), ylim = range(leo@coords[,2]), axes = TRUE)
points(leo@coords[,1], leo@coords[,2], col = "red", pch = ".")


#select summer breeding area
leo = leo[leo@coords[,1] < -106 & leo@coords[,2] > 53,]

# project leo to UTM, but keep leoGeo for lat-long
leoGeo = leo
leo = spTransform(leo, CRS = "+proj=utm +zone=13 +datum=WGS84")

# The recurse package supports both move objects and data frames
# here we demonstrate creating a data frame with the x, y, time, and id
leo.df = as.data.frame(leo@coords)
leo.df$time = leo$timestamp 
leo.df$id = "leo" 

# information to select radius
# calculate time and distance of each step
timeDiff = diff(as.numeric(leo.df$time) / 3600) # in hours
distDiff = sqrt(diff(leo.df$coords.x1)^2 + diff(leo.df$coords.x2)^2)

# step lengths eliminating those that cross the year boundary (~ over 6 months apart)
summary(distDiff[timeDiff < 4000]) 

# gaps (over 4000 are those that cross year boundary)
table(timeDiff)
sum(timeDiff > 1) / length(timeDiff)

#-------------------------------------
# Caluculate revisits in 50m radius
#-------------------------------------

# 50 m radius to find nesting/roosting sites (be careful of GPS error, approx. +/- 18m)
leovisit50 = getRecursions(leo.df, 50)  

# plot results
par(mfrow = c(2, 2), las = 1, bty = "l")
plot(leovisit50, leo.df, legendPos = c(280000, 5990000)) 
plot(leovisit50, leo.df, xlim = c(315000, 345000), ylim = c(5930000, 5960000))

hist(leovisit50$revisits, xlab = "Revisitations", main = "", col = "gray70")
hist(leovisit50$revisitStats$timeInside, xlab = "Time inside (h)", main = "", col = "gray70")

par(mfrow = c(1,1))

# where are the most frequently revisited locations?
revisitThreshold = 75

leoGeo.map.df = as(leoGeo,'data.frame')
leoGeo.map.df$revisits = leovisit50$revisits

# set zoom = 8 to see entire breeding area or zoom = 13 to see frequently visited places
map.leoGeo = qmap(bbox(extent(leoGeo[leovisit50$revisits > revisitThreshold,])), 
				  zoom = 11, maptype = "road", legend = "topright")

print(map.leoGeo + 
	  	geom_point(data = leoGeo.map.df[order(leoGeo.map.df$revisits),], 
	  			   aes(x = coords.x1, y = coords.x2, color = revisits), 
	  			   size = 0.75, alpha = 0.5)
	  + scale_color_viridis(name = "revisits", option = "magma", trans = "log10", 
	                        breaks = c(1, 2, 10, 50, 200))
)

#-------------------------------------
# environmental covariate
#-------------------------------------

# we use distance to roads as an environmental covariate to look at revists with

# check that road distance data frame is in same order as leo
all(leo$event.id == roadDist$event.id)

# summary
summary(roadDist$distToRoad)
sum(roadDist$distToRoad < 1000) / length(roadDist$distToRoad)
max(roadDist$distToRoad[leovisit50$revisits > revisitThreshold])
max(leovisit50$revisits[roadDist$distToRoad > 1000])

# add revists and lat/long to road distance data frame 
roadDist$revisits = leovisit50$revisits
roadDist$latitude = leoGeo@coords[,1]
roadDist$longitude = leoGeo@coords[,2]

print(map.leoGeo + 
        geom_point(data = roadDist, 
                   aes(x = latitude, y = longitude, color = distToRoad), 
                   size = 0.75, alpha = 0.8)
      + scale_color_viridis(name = "dist to road (m)")
)

roadvists = ggplot(roadDist, aes(x = distToRoad, y = revisits)) +
  xlab("distance to road (m)") + ylab("revisits") +
  geom_point(alpha = 0.2, col = "black") +
  theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
    axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid') )
roadvists


#-------------------------------------
# clusters
#-------------------------------------

# assign the most frequently visited locations to 5 clusters
leo.df.subset = leo.df[leovisit50$revisits > revisitThreshold,]
cluster = fanny(leo.df.subset[,c("coords.x1", "coords.x2")], k = 5)
leo.df.subset$cluster = as.factor(cluster$clustering)


# get centroid of each cluster, and also convert back to lat-long to plot on map
clusters = data.frame(x = tapply(leo.df.subset[,"coords.x1"], cluster$clustering, mean), 
					  y = tapply(leo.df.subset[,"coords.x2"], cluster$clustering, mean))

clusters = cbind(clusters, as.data.frame(spTransform(
	SpatialPoints(coords = clusters, proj4string = CRS("+proj=utm +zone=13 +datum=WGS84")), 
	CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))))
names(clusters)[3:4] = c("lon", "lat")

map.leoGeo = qmap(bbox(extent(leoGeo[leovisit50$revisits > revisitThreshold,])), 
				  zoom = 13, maptype = "road", legend = "topright")

print(map.leoGeo + 
	  	geom_point(data = leoGeo.map.df, aes(x = coords.x1, y = coords.x2), 
	  			   size = 0.5, alpha = 0.2, col = "gray")
	  + geom_point(data = subset(leoGeo.map.df, leovisit50$revisits > revisitThreshold), 
	  			 aes(x = coords.x1, y = coords.x2), 
	  			 size = 0.5, alpha = 0.3, color = "red")
	  + geom_point(data = clusters, aes(x = lon, y = lat),
	  			 size = 3, shape = 49:53)
) 


#-------------------------------------
# re-run revisits for just 5 clusters 
#-------------------------------------

# assign year as id to treat years separately
# again use 50m radius
leo.df$id = year(leo.df$time)
leositevisit = getRecursionsAtLocations(leo.df, clusters[,1:2], 
										radius = 50, threshold = 0.1) 

leostats = leositevisit$revisitStats

# to exclude very short duration visits
leostats = leostats[leostats$timeInside > 0.1,]


# convert to local time (timestamp is UTC time) and calculate some convience variables
leostats$local.entranceTime = with_tz(leostats$entranceTime, tz = "Canada/Saskatchewan")
leostats$local.exitTime = with_tz(leostats$exitTime, tz = "Canada/Saskatchewan")
leostats$year = year(leostats$entranceTime)
leostats$doy_enter = yday(leostats$local.entranceTime)
leostats$hour_enter = hour(leostats$local.entranceTime)
leostats$time_enter = hour(leostats$local.entranceTime) + minute(leostats$local.entranceTime) / 60 + second(leostats$local.entranceTime) / 60 / 60
leostats$time_exit = hour(leostats$local.exitTime) + minute(leostats$local.exitTime) / 60 + second(leostats$local.exitTime) / 60 / 60
leostats$overnight = factor(as.logical(yday(leostats$local.exitTime) - yday(leostats$local.entranceTime)))
levels(leostats$overnight) = c("no", "yes")
leostats$site = paste("site", leostats$coordIdx) # for plotting


#-------------------------------------
# entrance/exit time of day 
#-------------------------------------

# calculate time of day of entrace and exit (important to use local time not UTC)
geoLocs = as.matrix(clusters[leostats$coordIdx, c("lon", "lat")])
leostats$light_enter = calculateTimeOfDay(geoLocs, leostats$local.entranceTime)
leostats$light_exit = calculateTimeOfDay(geoLocs, leostats$local.exitTime)

table(leostats$light_enter)
table(leostats$light_exit)

#-------------------------------------
# inter- and intra-annual visit patterns across 5 sites 
#-------------------------------------

# plot entrace time by year and day of year
ttime = ggplot(leostats, aes(x = doy_enter)) + 
	geom_histogram(binwidth = diff(range(leostats$doy_enter))/7, color = "darkgray", fill = "gray") +
	facet_grid(site ~ year) + 
	xlab("visit day of year") + ylab("revisit frequency") + 
	theme_classic()
ttime 


#-------------------------------------
# visit duration by entrance time of day
#-------------------------------------

binhour = ggplot(leostats, aes(x = time_enter, y = timeInside)) +
	geom_density2d(color = "black") + ylim(0, 24) + 
	scale_color_brewer(palette = "Dark2", name = 'overnight') +
	xlab("visit entrance time") + ylab("visit duration (h)") +
	geom_point(alpha = 0.2, aes(col = overnight)) +
	theme_classic() + theme(
		axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
		axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
		legend.justification = c(0, 1), legend.position = c(0, 1))
binhour

#-------------------------------------
# does visit duration vary through breeding season?
#-------------------------------------
bindoy = ggplot(leostats, aes(x = doy_enter, y = timeInside)) +
  geom_density2d(color = "black") + ylim(0, 24) + 
  scale_color_brewer(palette = "Dark2", name = 'overnight') +
  xlab("visit entrance day of year") + ylab("visit duration (h)") +
  geom_point(alpha = 0.2, aes(col = overnight)) +
  theme_classic() + theme(
    axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
    axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
    legend.justification = c(1, 1), legend.position = c(1,1)) +
  facet_grid(~ site)
bindoy

#-------------------------------------
# residence time
#-------------------------------------

# calculate approximate residence time by rounding to nearest hour (data is hourly)
begin = round(leostats$time_enter)
duration = round(leostats$timeInside)
hours = data.frame(site1 = rep(0, 24),
				   site2 = rep(0, 24),
				   site3 = rep(0, 24),
				   site4 = rep(0, 24),
				   site5 = rep(0, 24))

for (i in 1:nrow(leostats))
{
	#idxs will be 0-23, so add 1 for array indexing
	idxs = (begin[i] + 0:(duration[i])) %% 24
	hours[idxs + 1, leostats$coordIdx[i]] = hours[idxs + 1, leostats$coordIdx[i]] + 1
}

# reformat data for plotting with reshape2
hours$hour = factor(row.names(hours), levels = as.character(1:24), ordered = TRUE)
hours.melt = melt(hours, value.name="Count", variable.name="site", na.rm = TRUE)

restime = ggplot(hours.melt, aes(x = hour, y = Count)) +
	geom_bar(stat = "identity", color = "darkgray", fill = "gray") +
	xlab("hour of day") + ylab("total hours") + 
	scale_x_discrete(breaks = 1:24, 
					 labels = c("1", "", "", "", "", "6", "", "", "", "", "", "12", "", "", "", "", "", "18", "", "", "", "", "", "24")) +
	theme_classic()
restime + facet_grid(~ site) 

#-------------------------------------
# time since last visit
#-------------------------------------

# seems to be a difference between sites 1-3 and 4-5 for residence time
leostats$siteGroup = ifelse(leostats$coordIdx < 4, "sites 1-3", "sites 4-5")

# does time spend vary based on recency of last visit? i.e. less than a day, two days, or longer
leostats$timesince = cut(leostats$timeSinceLastVisit, 
						 breaks = c(0, 24, 48, 2000), labels = c("<24h", "24-48h", ">48h"))

boxtimesince = ggplot(na.omit(leostats), aes(x = timesince, y = timeInside, fill = siteGroup)) + 
	geom_boxplot() + ylim(0, 50) + 
	scale_fill_manual(name = "", values = c("white", "gray")) +
	xlab("time since last visit") + ylab("visit duration (h)") +
	theme_classic() +
	theme(legend.justification = c(1, 1), legend.position = c(0.9, 0.9) ) 
boxtimesince 
