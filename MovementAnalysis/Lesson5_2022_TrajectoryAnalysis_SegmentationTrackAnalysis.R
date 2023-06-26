#########################################
###           AniMove 2022            ###    
### Script by Kami Safi & Anne Scharf ###
#########################################

library(maptools)
library(move)
library(lubridate)
library(mapdata)
library(scales)
setwd("/home/ascharf/Documents/Animove22/MovementAnalysis/data/")


#########################
## BURST A TRAJECTORY ###
#########################
# we define night as the time when the sun at the beginning and end of a segment is 6 degrees below the horizon
Leroy <- move(system.file("extdata","leroy.csv.gz",package="move"))
DayNight <- rep("Day", n.locs(Leroy)-1)
DayNight[solarpos(Leroy[-n.locs(Leroy)], timestamps(Leroy)[-n.locs(Leroy)])[,2] < -6 & 
           solarpos(Leroy[-1], timestamps(Leroy)[-1])[,2] < -6] <- "Night"
#assigning to each segment if it is during daytime or night
Leroy.burst <- move::burst(x=Leroy, f=DayNight)
Leroy.burst
str(Leroy.burst)

plot(Leroy.burst, type="l", col=c("red", "black"))
legend("bottomleft",legend=c("day","night"), col=c("red", "black"), lty=1)


#####################
###### CORRIDORS ####
#####################
# identify corridor behaviour, i.e. parallel fast movement. For details see ?corridor
LeroyCorr <- corridor(Leroy)
plot(LeroyCorr, type="l", xlab="Longitude", ylab="Latitude", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")


##############################
## Net Square Displacement ###
##############################
# load the move object of Sierit the stork
load("Sierit.RData")
# thinning the data to one position a day just to make the example code run faster. This is a "quick and dirty" way to thin a track, see also amt::track_resample
Sierit1day <- Sierit[!duplicated(round_date(Sierit@timestamps,"24 hour"))]
# first inspect the data
par(mar=c(0,0,0,0))
ext <- bbox(Sierit1day)+c(-5,-5,5,5)
map('worldHires', xlim = ext[1, ], ylim = ext[2, ])
lines(Sierit1day, col="red")
points(Sierit1day, col=alpha("red",0.3), pch=20)

library(adehabitatLT)
## for migrating birds for which we have the nest location, NSD often gives a nice overview of what is happening
SieritNSD <- (spDistsN1(coordinates(Sierit1day), coordinates(Sierit1day)[1,],longlat=T))^2
plot(Sierit1day$timestamp, SieritNSD, type="l",xlab="Time", ylab="Net square distance (Km²)", main="Sierit")
points(Sierit1day$timestamp, SieritNSD, pch=20, cex=0.7)

# finding break points in the NSD using the lavielle function
lvNSD <- lavielle(SieritNSD, Lmin=9, Kmax=12, type = "var") 
brksNSDl <-  findpath(lvNSD, 12, plotit = F)
brksNSD <- unlist(brksNSDl)
plot(Sierit1day$timestamp, SieritNSD, type="l", xlab="", ylab="NSD (Km²)", main="NSD")
points(Sierit1day$timestamp, SieritNSD,pch=20,cex=.4)
abline(v = Sierit1day$timestamp[brksNSD], col = "red")

ts_cutNSD <- Sierit1day$timestamp[brksNSD]
# assigning unique ID to each segment
for(x in seq(1,length(ts_cutNSD),2)){
  Sierit1day$nsdcat[Sierit1day@timestamps>=ts_cutNSD[x] & Sierit1day@timestamps<=ts_cutNSD[x+1]] <- paste0("Seg_",x)
}
# bursting move object by segment ID
SieritNSD <- move::burst(Sierit1day,f=Sierit1day$nsdcat[-1])
plot(SieritNSD, type="l", col=primary.colors(length(levels(SieritNSD@burstId)),steps=2), main="NSD")
points(SieritNSD,pch=20, col=primary.colors(length(levels(SieritNSD@burstId)),steps=2)[SieritNSD@burstId])

## (end of segmentation section) ##


#####################
## TRACK ANALYSIS ###
#####################
# calculate mean speed, distance, time, to do comparisons between individuals, sex, etc
library(move)
bats <- move("Parti-colored bat Safi Switzerland.csv")
medianSpeed<-unlist(lapply(speed(bats), median))
timeTracked<-unlist(lapply(timeLag(bats, units='days'), sum))
distanceTracked<-unlist(lapply(distance(bats), sum))
indData<-data.frame(medianSpeed, timeTracked, distanceTracked)
head(indData, 4)


## test for differences in distance traveled per day between sexes
# if the data are downloaded with the "getMovebankData" function the information on sex (if available) is included in the @idData
# the reference table can be now also downloaded via R with "getMovebankReferenceTable"
bats_ref <- read.csv("Parti-colored bat Safi Switzerland-reference-data.csv", sep=",", as.is=TRUE) 
bats_ref <- bats_ref[!is.na(bats_ref$animal.id),c("animal.id", "animal.sex")]
names(bats_ref) <- c("id", "sex")
bats_ref$id <- paste0("X", bats_ref$id)
indData$id <- row.names(indData)
indData <- merge(indData, bats_ref, by="id")
boxplot(log10(I(distanceTracked/timeTracked))~sex, data=indData, names=c("Females", "Males"), 
        ylab=expression(paste(Log_10, " of cumulative distance in m per day", sep="")))
# Note to plot: there may be a small difference 
t.test(log(I((distanceTracked/1000))/timeTracked)~sex,data=indData)
# but this difference does not seem to be significant
mod <- glm(sqrt(distanceTracked)~as.factor(sex)+timeTracked,data=indData)
par(mfrow=c(2,2))
plot(mod, ask=F)
summary(mod)
# model confirms


## test for differences in speed between sexes 
library(MASS)
boxplot(log(medianSpeed)~sex, data=indData, names=c("Females", "Males"), ylab="Median speed in m/s")
# Note to plot: males seem to travel faster than females
wilcox.test(medianSpeed~sex, data=indData)
# this time the difference seems to be significant

# boxcox is a nice function to find the right transformation of the data
bc <- boxcox(medianSpeed~as.factor(sex)+timeTracked, data=indData[-15,])
modII <- glm(I(medianSpeed^bc$x[which.max(bc$y)])~ as.factor(sex)+timeTracked, data=indData[c(-15),])
# plot(modII)
summary(modII)
# there is a difference between sexes in speed
# but!! this difference could be due to sampling frequency!


## changes in speed due to differences in sampling ##
library(scales)
library(adehabitatLT)
set.seed(7478)
# creating one track
r.track <- as(simm.crw(1:1000, 1, 0.99), "Move")
r.trackThin <- r.track[seq(1,nrow(coordinates(r.track)),3),]
# compare speed of the track, and speed of the same track but only taking every 3rd position
# non-parametric test to compate distributions
ks.test(sqrt(speed(r.track)), sqrt(speed(r.trackThin)))

hist(sqrt(speed(r.trackThin)), freq=F, xlim=c(0,2), breaks="FD", col=alpha("grey", 0.5), xlab="Square root transformed speed", main=NA)
hist(sqrt(speed(r.track)), freq=F, add=T, breaks="FD", col=alpha("white", 0.5))
legend("topleft", c("all positions", "every 3rd position"), fill=c("white", "grey"), bty="n")

# non-parametric test to compare means
wilcox.test(sqrt(speed(r.track)), sqrt(speed(r.trackThin)))
# parametric test to compate means
t.test(sqrt(speed(r.track)), sqrt(speed(r.trackThin)))



## using generalized linear mixed models
# - can include complex correlation structures
# - can include random factors
library(mgcv)
testDat <- data.frame(v=unlist(speed(bats)), 
                      id=trackId(bats)[!diff(c(0,as.numeric(trackId(bats))))],
                      time=(unlist(timestamps(bats))[!diff(c(0,as.numeric(trackId(bats))))])-(unlist(timeLag(bats, units="secs"))/2))
testDat <- testDat[testDat$v<20 & testDat$v>2,]
testDat <- merge(testDat, bats_ref, by="id")
cs1 <- corAR1(0.5, form = ~ time|id) # correlation structure

mod <- gamm(log(v)~as.factor(sex), random=list(id=~1), correlation=cs1, data=testDat, family=gaussian)
summary(mod$gam)
acf(residuals(mod$gam))
par(mfrow=c(2,2))
gam.check(mod$gam)



