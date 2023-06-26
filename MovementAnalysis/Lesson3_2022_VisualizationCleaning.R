#########################################
###           AniMove 2022            ###    
### Script by Kami Safi & Anne Scharf ###
#########################################

library(move)
setwd("/home/ascharf/Documents/Animove22/MovementAnalysis/data/")

###################################
#### MAPPING MOVEMENT DATA ########
###################################
bats <- move("Parti-colored bat Safi Switzerland.csv")

### basic plots ###
plot(bats)
plot(bats, xlab="Longitude", ylab="Latitude",type="b", pch=16, cex=0.5)

### plot on the world ###
library(mapdata)
library(scales)
map('worldHires', col="grey", fill=T)
points(t(colMeans(coordinates(bats))), col=alpha('red',0.5), pch=16)
points(t(colMeans(coordinates(bats))), col='cyan')

### plot on the world, zoomed in ####
(e<-bbox(extent(bats)*5))
# note here that the brackets around the assignment ensure that the result is also printed to the console
map('worldHires', xlim = e[1, ], ylim = e[2, ])
points(bats)
lines(bats)

### plot on google background ####
library("ggmap")
library("mapproj")
# coerce move object to a data frame
bats_df <- as.data.frame(bats)
# request map data from google
m <- get_map(e, zoom=9, source="google", maptype="terrain")
# plot the map and add the locations separated by individual id
ggmap(m)+geom_path(data=bats_df, aes(x=location.long, y=location.lat, colour=trackId))
# we can also add a scalebar
library(ggsn)
xylim <- as.numeric(attributes(m)$bb)
ggmap(m)+geom_path(data=bats_df, aes(x=location.long, y=location.lat, colour=trackId))+
  scalebar(x.min = xylim[2], x.max = xylim[4],
           y.min = xylim[1], y.max = xylim[3],
           dist = 50, dist_unit="km", transform=T, model = 'WGS84',anchor=c(x=11.6,y=45.6),st.size=3)


### interactive maps ###
library(mapview)
# with a movestack just to inspect data
mapview(bats)

# transform to e.g. spdf and sldf for more options
batsSPDF <- as.data.frame(bats)
coordinates(batsSPDF) <- ~location.long + location.lat
projection(batsSPDF) <- projection(bats)

batsLsL <- lapply(split(bats), function(x){Lines(list(Line(coordinates(x))),ID=namesIndiv(x))})
batsSLDF <- SpatialLinesDataFrame(SpatialLines(batsLsL,proj4string = bats@proj4string), data=idData(bats) ,match.ID = F)

mapview(batsSPDF,zcol="individual.local.identifier")+mapview(batsSLDF,zcol="individual.local.identifier", legend=F)

library(tmap)
# with a movestack just to inspect data
tmap_mode("view")
tm_shape(bats)+tm_dots()

# transform to e.g. sf class for more options
library(sf)
bats_SFp <- bats_df%>%st_as_sf(coords = c("location.long", "location.lat"), crs = crs(bats))%>%st_cast("POINT")
tmap_mode("view")
tm_shape(bats_SFp)+tm_dots(col = "individual.local.identifier")

## see this book "Making Maps with R" for more details on making maps in R: 
## https://bookdown.org/nicohahn/making_maps_with_r5/docs/introduction.html


####################################################
####### REMOVING OUTLIERS BASED ON MAPPING #########
####################################################

load("buffalo_cleaned.Rdata") # buffalo

## Create gray scale
buffaloGray<-gray((nrow(idData(buffalo))-1):0/nrow(idData(buffalo)))
## Plot with gray scale
plot(buffalo, col=buffaloGray, xlab="Longitude", ylab="Latitude")

## get the position of the coordinate that has the max longitude
which.max(coordinates(buffalo)[,1])
## drop the point with the largest coordinate values
buffalo <- buffalo[-which.max(coordinates(buffalo)[,1])]
plot(buffalo, col=buffaloGray, xlab="Longitude", ylab="Latitude")
## save the clean dataset for the following days
save(buffalo, file="buffalos.Rdata")


#################################################
### TEMPORAL ORGANIZATION OF THE TRAJECTORIES ###
#################################################
## number of locations
n.locs(bats)

## time lag between locations
timeLags <- timeLag(bats, units='hours') # important: always state the units!

## distribution of timelags
timeLagsVec <- unlist(timeLags)
summary(timeLagsVec)
hist(timeLagsVec, breaks=50, main=NA, xlab="Time lag in hours")
arrows(24.5,587.5,20.7,189.7, length=0.1)
arrows(49.5,587.5,45.7,189.7, length=0.1)

## distribution of timelags shorter than 1h
hist(timeLagsVec[timeLagsVec<1], breaks="FD", main=NA, xlab="Time lag in hours")
## count of locations per timebin
summary(as.factor(round(timeLagsVec, 4)), maxsum=5)

## nb locations per hour
ts <- timestamps(bats)
library('lubridate')
#transform timestamps into local time of study for better interpretation
tsLocal <- with_tz(ts, tzone="Europe/Zurich")
tapply(tsLocal, hour(tsLocal), length)

## nb locations per month and hour
tapply(tsLocal, list(month(tsLocal),hour(tsLocal)), length)


#########################################
### SPATIAL ORGANIZATION OF THE TRACK ###
#########################################

### distance between locations ###
dist <- unlist(distance(bats))
summary(dist)
hist(dist)

### speed between locations ###
speeds <- unlist(speed(bats))
summary(speeds)
hist(speeds, breaks="FD")

## have a look at the realistic speeds (e.g.<20m/s)
speedsRealistic <- speeds[speeds<20]
speedVsTimeLag <- data.frame(timeLag=timeLagsVec, speeds=speeds)
speedVsTimeLag <- speedVsTimeLag[speedVsTimeLag$timeLag<10 & speedVsTimeLag$speeds<20,]
plot(speedVsTimeLag$timeLag, speedVsTimeLag$speeds, xlab='Time lag (hours)', ylab='Speed (m/s)', pch=19)
hist(speedsRealistic, main=NA, xlab="Speed in m/s", breaks="FD", ylab="Frequency")
## this is the common shape for speeds


## identify the segments with the high and low speeds ##
## => start <= ##
library(classInt)
# select bat X191
bat191 <- bats[["X191"]]
# store speed
v <- speed(bat191)
# find 49 breaks in the speed vector
my.class <- classIntervals(unlist(v),n=49,style="equal")
# assign colours in 50 shades of grey
my.pal <- findColours(my.class,grey(0:49/49))

# make a data frame that is later used to draw the segments
# with x0, y0, x1, y1 coordinates
segdat <- as.data.frame(cbind(coordinates(bat191)[-n.locs(bat191),], coordinates(bat191)[-1,]))
# add colours to the data frame
segdat$col <- as.vector(my.pal)
# add the speed
segdat$v <- v
# sort the data frame by increasing speed
# this will make sure high speed segments are plotted on top
# of low speed segments
segdat <- segdat[order(segdat$v),]

# change the margins of the plot region
par(mar=c(5,4,4,5))
# create a plot with the appropriate size but no points
plot(bat191, xlab="Longitude", ylab="Latitude", main="Speed in m/s", asp=1, type="n")
# get the size of the plot region
u <- par("usr")
# draw a rectangle of the size of the plot region and colour it 
# with the colour that corresponds to the median speed colour
rect(u[1], u[3], u[2], u[4], col = segdat$col[length(v[v <= median(v)])])
# now draw the segments
segments(segdat[,1],segdat[,2], segdat[,3], segdat[,4], col=segdat[,5], lwd=3)
# add a legend
plot(raster(nrows=1, ncols=length(v), vals=v), legend.only=TRUE, col=grey(0:49/49), legend.mar=4.5)
## => end <= ##


### direction of movement / azimuth / heading of movement ###
## NOTE: heading or bearing are mostly refer to the direction of body axis
direction <- angle(bats) ## Angles in degrees relative to the N
summary(unlist(direction))
hist(unlist(direction),  breaks=18, xlab="Direction of movement", main=NA)
# they seem to go into all directions

### turning angles ###
turnAngles <- turnAngleGc(bats) ## Angles in degrees relative to the previous step
hist(unlist(turnAngles), breaks=18, xlab="Turning Angle", main=NA)
## this shape can indicate movement along a linear structure

turnAnglesBuf <- unlist(turnAngleGc(buffalo))
hist(turnAnglesBuf)
## this is the common shape for turning angles


############################
######## MISSED FIXES ######
############################

library("lubridate")
# load the data
leroy <- move(system.file("extdata","leroy.csv.gz",package="move"))
# get the time stamps of entries without location
pattrn <- data.frame(time=leroy@dataUnUsedRecords$timestamp, status='Not Successful')
# add the time stamps of positions obtained
pattrn <- rbind(pattrn, data.frame(time=timestamps(leroy), status='Successful'))
# change time to local time zone
pattrn$time<-with_tz(pattrn$time, tz="America/New_York")

# Load ggplot library for plotting
library("ggplot2")
# Plot histogram that is filled out for proportions and is binned per hours
ggplot(pattrn, aes(x=hour(time), fill=status))+
  geom_histogram(binwidth=1, position='fill')+scale_fill_grey()

## histogram filled out for proportions binned per day
ggplot(pattrn, aes(x=time, fill=status))+
  geom_histogram(binwidth=24*60*60, position='fill')+scale_fill_grey()



                                     