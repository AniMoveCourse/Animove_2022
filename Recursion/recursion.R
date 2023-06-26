#########
# Animove recursion exercise
# Chloe Bracis and Thomas Mueller
# chloe.bracis@gmail.com

#load libraries
library(recurse)
library(move)
library(ggmap)
library(RgoogleMaps)
library(lubridate)
library(raster)
library(scales)


############################
#two function for color palettes

#get continuous palette from blue (low values) to red (high values)
getContinuousPalette = function(n, alpha = 1)
{
  cols = alpha(brewer_pal(palette = "RdYlBu", direction = -1)(9), alpha)
  return( gradient_n_pal(cols)(seq(0, 1, length = n)) )
}

######################
#discrete sequential palette, e.g. for one color per day 
getDiscreteSequentialPalette = function(n, alpha = 1)
{
  return(alpha(rainbow(n), alpha))
}
####################################

############################
# read movebank elephant "Habiba" from Wall et al. 2014


elephants =  split(move("Elliptical Time-Density Model (Wall et al. 2014) African Elephant Dataset (Source-Save the Elephants).csv"))
habiba = elephants$Habiba
plot(habiba)



#summarize data
summary(habiba)
range(habiba$study.local.timestamp)
diff(range(habiba$utm.easting))
diff(range(habiba$utm.northing))
##################################


#Color code data day by day
uniqueDays = unique(day(habiba$timestamp))
dayIdx = sapply(day(habiba$timestamp), function(x) which(x == uniqueDays))
move::plot(habiba, 
           col = getDiscreteSequentialPalette(length(uniqueDays))[dayIdx])
# -> it appears that the looping behavior occurs at night
#####################################



#####################################
#Examine the Revistation data
######################

# project to a equidistant projection (move's spTransform uses aeqd by default)
# don't want to do this with latlong because length of a degree changes with latitude!
habibaAEQD = spTransform(habiba, center = TRUE)

# calculate recursions for a 500m Radius
Radius = 500# 500m 

habibavisits = getRecursions(habibaAEQD , radius = Radius) 
names(habibavisits)
#habibavisits = getRecursions(habiba.df , radius = Radius) 

# examine output object
str(habibavisits)
head(habibavisits$revisits)

#plot revisits
hist(habibavisits$revisits, breaks = 10, col = "darkgrey", border = NA, main = "", xlab = "Revisits (radius 500m)")

#plot revisits on trajectory
plot(habibavisits, habibaAEQD, alpha = 1, legendPos = c(-3000, 4000))
# draw circle as a reference for scale
drawCircle(max(habibaAEQD@coords[,1]) - Radius, min(habibaAEQD@coords[,2]) + Radius, radius = Radius)
######################################


######################################
# plot in google maps
# note: google maps expects latlong coordinates
habiba.map.df = as(habiba,'data.frame')
map.habiba <- qmap(bbox(extent(habiba)), zoom = 13, maptype = 'hybrid', legend="topright")
print(map.habiba + 
        geom_path(data = habiba.map.df, aes(x = coords.x1, y = coords.x2), color = "white", size = 0.3) + 
        geom_point(data = habiba.map.df, aes(x = coords.x1, y = coords.x2), 
                   color = getContinuousPalette(max(habibavisits$revisits), 0.5)[habibavisits$revisits]))



###############################################################
# there is more than just the revisits in the recursion object:
# accessing the data frame for the recursion statistics
habibavisits$revisitStats[1:15,]

#################################
# examine first passage time
# the first visits passes through the center of the circle, thus representing a first passage time
# the first visit at each location equals the first passage time
# this only works in teh recurse package if you use single individuals, not available for multiple individuals
habibavisits$firstPassageTime = habibavisits$revisitStats$timeInside[habibavisits$revisitStats$visitIdx==1]


hist(as.numeric(habibavisits$firstPassageTime), breaks = 20, col = "darkgrey", border = NA, main = "", xlab = "First passage (hrs)")
# the histogram indicates a bimodal distribution potentially indicating two different behaviors where


#plot locations with first passage larger or smaller 6 hrs on map
cutOff = 6
print(map.habiba + 
        geom_path(data = habiba.map.df, aes(x = coords.x1, y = coords.x2),color = "white", size = 0.3) + 
        geom_point(data = habiba.map.df, aes(x = coords.x1, y = coords.x2), 
                   color=alpha(ifelse(habibavisits$firstPassageTime >cutOff, "blue", "grey"),.5)))


#it seems there is a different behavior at sections of the loop, maybe related to daytime if loops are daily
boxplot(as.numeric(habibavisits$firstPassageTime) ~ hour(habiba$timestamp), 
		outline = FALSE, col = "grey", xlab = "Daytime (hrs)", ylab = "First passage (hrs)")



########################
# examine residence/utilization time, the sum of all visits around a focal point
head(habibavisits$residenceTime)


hist(as.numeric(habibavisits$residenceTime), breaks = 20, col = "darkgrey", border = NA, main = "", xlab= "Utilization time (hrs)")
# there seems to be a bimodal distribution, separated at about 20 hrs total visit time
colors = as.character(cut(habibavisits$residenceTime, breaks = c(0, 10, 20, 30, 40), labels = c("grey", "gold", "darkorange", "firebrick3")))
print(map.habiba + 
        geom_path(data = habiba.map.df, aes(x = coords.x1, y = coords.x2), color = "white", size = 0.3) + 
        geom_point(data = habiba.map.df, aes(x = coords.x1, y = coords.x2), 
                   color = alpha(colors, 0.5)))

#it seems there is a different behavior with loops, maybe related to nighttime if loops are daily
boxplot(as.numeric(habibavisits$residenceTime) ~ hour(habiba$timestamp), 
        outline = FALSE, col = "grey", xlab = "Daytime (hrs)", ylab = "Residence passage (hrs)")


########################
#look at NDVI at these areas to explore relation to resources
# note: the NDVI raster is in UTM coordinates
ndvi = raster("MOD13Q1.A2014033.250m_16_days_NDVI.tif") / 10000
habiba.ndvi = data.frame(habiba[,c("utm.easting", "utm.northing", "timestamp")] )[,1:3] # ignore extra cols move package adds
habiba.ndvi$ndvi = extract(ndvi, habiba.ndvi[,c("utm.easting", "utm.northing")])


# plot NDVI with data
plot(ndvi)
lines(habiba.ndvi$utm.easting, habiba.ndvi$utm.northing, col = "white")
points(habiba.ndvi$utm.easting, habiba.ndvi$utm.northing, 
       col = alpha(colors, 0.3), pch = 16, cex = 0.7)
drawCircle(max(habiba.ndvi$utm.easting) - Radius, min(habiba.ndvi$utm.northing) + Radius, radius = Radius)

# check if there is any difference in NDVI between residence time categories
cutOff = 20
graphics::boxplot(habiba.ndvi$ndvi ~ habibavisits$residenceTime > cutOff, 
				  outline = FALSE, col = "grey", notch = FALSE, 
				  xlab = paste0("Utilization > ", cutOff," hrs"), ylab = "NDVI")



# check if there is any difference in NDVI between times of day
graphics::boxplot(habiba.ndvi$ndvi ~ hour(habiba$timestamp), 
                  outline = FALSE, col = "grey", xlab = "Time of day (hrs)", ylab = "NDVI")



# look at revisitation
hist(as.numeric(habibavisits$revisitStats$timeSinceLastVisit / 24), freq = TRUE, 
     xlab = "Time since last visit (days)" , col = "darkgrey", border = NA, main = "")

returnsAfterOneWeek = as.vector(na.omit(habibavisits$revisitStats$coordIdx[as.numeric(habibavisits$revisitStats$timeSinceLastVisit / 24) > 7]))

print(map.habiba + 
        geom_path(data = habiba.map.df, aes(x = coords.x1, y = coords.x2), color = "white", size = 0.3) + 
        geom_point(data = habiba.map.df[returnsAfterOneWeek, ], aes(x = coords.x1, y = coords.x2), 
                   color = alpha("blue", 0.2)))


# look at revisitation times at shorter time intervals
hist(as.numeric(habibavisits$revisitStats$timeSinceLastVisit / 24), freq = TRUE, 
   xlab = "Time since last visit (days)" , xlim = c(0, 3), ylim = c(0, 400), 
   breaks = 70, col = "darkgrey", border = NA, main = "")
# there seems to be a hint for periodivity - possible suggests periodogram analyses


# what if we pick a smaller radius?
Radius = 100 # 100m 

habibavisits = getRecursions(habibaAEQD, radius = Radius) 

hist(habibavisits$revisits, breaks = 10, col = "darkgrey", border = NA, main = "", xlab = "Revisits (50m radius)")

plot(habibavisits, habibaAEQD, alpha = 1, legendPos = c(-3000, 4000))
drawCircle(max(habibaAEQD@coords[,1]) - Radius, min(habibaAEQD@coords[,2]) + Radius, radius = Radius)

# check out the package vignette for an example of testing different radii


















