#########################################
###           AniMove 2022            ###    
### Script by Kami Safi & Anne Scharf ###
#########################################

library(move)
library(lubridate)
library(circular)
library(fields)
library(mapdata)
library(scales)
setwd("/home/ascharf/Documents/Animove22/MovementAnalysis/data/")

############################################
#### AZIMUTH, TURNING ANGLE and SPEED ######
############################################
## load and plot Leo
Leo <- move("Leo-65545.csv.gz")
map('worldHires', xlim = Leo@bbox[1, ]+c(-5,5), ylim = Leo@bbox[2, ]+c(-5,5),col="grey", fill=T)
lines(Leo, col="firebrick", lwd=1.5)

## get locations per year and month
tapply(timestamps(Leo), list(year(timestamps(Leo)), month(timestamps(Leo))), length)

## removing years 2012-13 because there is a large gap
Leo <- Leo[year(Leo$timestamp)<2012,]


###### azimuth distribution per season (circular) ##### 
# Define the categorization. Numbers refer to month
categories<- c("1"="Wintering", "2"="Wintering", 
               "3"="Wintering", "4"="North migration", 
               "5"="North migration", "6"="Breeding", 
               "7"="Breeding", "8"="Breeding", 
               "9"="South migration", "10"="South migration", 
               "11"="Wintering", "12"="Wintering")
# assign the categories to a new variable based on the timestamp
Leo$cat<-factor(categories[month(timestamps(Leo))],
                   levels=c("South migration", "Wintering", "North migration", "Breeding"))

# assign NA to segments where subsequent locations fall in different seasons
Leo$cat[c(Leo$cat[-n.locs(Leo)]!=Leo$cat[-1], TRUE)] <- NA

# store the information in a new data frame
azimuth <- data.frame(D=angle(Leo),
                      V=speed(Leo), 
                      Season=Leo$cat[-1])
# Define the direction as a circular
azimuth$Dcirc<-as.circular(azimuth$D, 
                 rotation="clock", 
                 units="degrees", 
                 type="angles", 
                 modulo="asis", 
                 zero=0, 
                 template="geographic")
# select segments above 2 m/s, we are only interested in segments when Leo is moving, and not the stationary error
azimuth <- azimuth[azimuth$V>=2,]
# remove missing values
azimuth <- azimuth[complete.cases(azimuth),]
# define a vector that is used to set the order of plotting
seasons <- levels(Leo$cat)
# change margins of plot
par(mar=rep(1,4))
# plot all the azimuths
plot(azimuth$Dcirc, stack=T, shrink=1.6, pch=16, sep=0.05, col="grey")
# loop through seasons and plot a line denisty per season
for(i in 1:length(seasons)){
  # subset the azimuth
  x <- azimuth[azimuth$Season==seasons[i],'Dcirc']
  # calculate density and plot as a line
  lines(density(x, bw=180, kernel="vonmises"), lwd=2, lty=i)
  # draw an arrow showing mean and resultant length
  arrows.circular(mean(x), y=rho.circular(x), lwd=2, length=0.1, lty=i)
}
# add a legend
legend("bottomleft", lty=c(1,2,3,4), seasons, bty="n", cex=0.85)
# Plot explanation (only contains directional information): 
# - grey histogram: number of locations
# - lines: circular density function
# - arrows: mean direction; arrow length: mean resultant length, a measure of concentration


## speed ~ azimuth scatter plot ###
library(scales)
plot(speed(Leo)~angle(Leo),
     ylab="Speed in m/s", xlab="Azimuth in degrees", type="n", bty="n")
points(speed(Leo)~angle(Leo), ylim=c(0,20),
     ylab="Speed in m/s", xlab="Azimuth in degrees", pch=16, col=alpha("black", 0.3))


## wind rose of azimuth and speed per season ###
par(list(mfrow=c(2,2)))
for(i in seasons){
  windrose(x=azimuth[azimuth$Season==i,'Dcirc'], 
           y=azimuth[azimuth$Season==i,'V'],
           main=i, plot.mids=T, cir.ind = 0.2, 
           mids.size=1, increment=5, bins=36, 
           fill.col=grey(seq(1,0, length.out=6)), 
           shrink=1)
}
par(mfrow=c(1,1))
# Note to plot:
# - migration has directionality with high speeds
# - breeding has only low speeds
# - wintering seems to contain some migration. But we just took month to categorize.


## wind rose of turning angle and speeds ###
turn <- data.frame(angle=turnAngleGc(Leo))
v <- speed(Leo)
turn$Vm <- rowMeans(cbind(as.numeric(v)[-1], as.numeric(v)[-length(as.numeric(v))]))
segSeason <- rowMeans(cbind(as.numeric(as.factor(Leo$cat))[-1], as.numeric(as.factor(Leo$cat))[-length(Leo$cat)]))
turn$season <- rowMeans(cbind(segSeason[-1], segSeason[-length(segSeason)]))
turn$season <- factor(turn$season, labels=levels(as.factor(Leo$cat)))
turn <- turn[turn$Vm>2,]
turn <- turn[!is.na(turn$season),]
# turning angle goes from -180 to 180, so we have to transform them into 0-360, to be able to plot them on a windrose
angle360 <- turn$angle[!is.na(turn$angle)]
angle360[angle360<0] <- angle360[angle360<0]+360
par(bty="n")
windrose(as.circular(angle360, 
                     rotation="clock", 
                     units="degrees", 
                     type="direction", 
                     modulo="asis", 
                     zero=0, 
                     template="none"), 
         y=turn$Vm[!is.na(turn$angle)],
         plot.mids=T, cir.ind = 0.2, mids.size=1,
         increment=5, bins=72, fill.col=grey(seq(1,0, length.out=6)),
         main="Turning angle and speed", tcl.text=-0.07)
# Note to plot:
# - high speeds only have low turning angles (most animals can't make high turns when going fast)
# - low speeds have all turning angles


## speed ~ turning angle per season ###
library(scales)
par(mfrow=c(2,2))
for(i in unique(turn$season)){
  plot(Vm~angle, data=turn[turn$season==i,], 
       xlab="Turning angle", ylab="Speed in m/s",
       bty="n", pch=16, ylim=c(0,max(turn$Vm)),
       col=alpha("black", 0.3), main=i)
}
par(mfrow=c(1,1))
# Note to plot:
# - during migration movement is directional and speeds are higher
# - in breeding speeds are low and movement is in all directions
# - here again we see that wintering contains some migration



#####################################################
## MOVEMENT PROCESS AND EFFECTS ON PATH METRICS #####
#####################################################
library(adehabitatLT)
library(move)
# simulate 500 tracks with 1000 steps, with different levels of correlation of turning angles (r) and convert into move object
sets <- sort(rep((0.8 + log10(c(seq(1,100, length.out=10)))/10)[1:9],500))
rCRW <- lapply(lapply(sets, simm.crw, date=1:1000, h=1), as, "Move")
# calculate NSD for all tracks, from origin of trajectory
rNSD <- unlist(lapply(lapply(lapply(rCRW, coordinates), spDistsN1, pt=c(0,0)), "^", 2))

mNSD <- tapply(rNSD, list(sort(rep(sets,1000)), rep(1:1000, length(sets))), mean)
par(mar=c(5, 4, 4, 4) + 0.1)
plot(0,0, type="n", xlim=c(0,1300), ylim=c(0, max(mNSD)),
     bty="n", xlab="Step", ylab="Net square distance", xaxt="n")
axis(1, at=c(0,200,400,600,800,1000))
test <- apply(as.matrix(mNSD), 1, lines, x=1:1000)
text(cbind(rep(c(1250, 1100), length.out=length(row.names(mNSD))), mNSD[,ncol(mNSD)]), 
     paste("r=", as.character(round(as.numeric(row.names(mNSD)),3)),sep=""), cex=0.5)
# Note to plot:
# - slope: how fast it gets away from the origin
# - r=0.99: very correlated, movement almost in a steady direction, it's walking away from the origin. For example:
plot(rCRW[[4500]], pch=19)
points(rCRW[[4500]][1], col="green3", pch=20)
points(rCRW[[4500]][1000], col="red", pch=20)
# - r=0.8: fairly correlated. movement is wiggly. For example:
plot(rCRW[[1]], pch=19)
points(rCRW[[1]][1], col="green3", pch=20)
points(rCRW[[1]][1000], col="red", pch=20)



#######################################
### NET SQUARE DISPLACEMENT (NSD) #####
#######################################
layout(matrix(c(1,1,2,3), ncol=2, byrow=T))
LeoNSD <- (spDistsN1(coordinates(Leo), coordinates(Leo)[1,],longlat=T))^2
plot(Leo$timestamp, LeoNSD, type="l",
     xlab="Time", ylab="Net square distance (Km²)", main="All data")

leoBreed08 <- Leo[which(Leo$cat=="Breeding" & year(Leo$timestamp)==2008),]
leoBreed08NSD <- (spDistsN1(coordinates(leoBreed08), coordinates(leoBreed08)[1,],longlat=T))^2
plot(leoBreed08$timestamp,leoBreed08NSD, type="l",
     xlab="Time", ylab="Net square distance (Km²)", main="Breeding 2008")

leoWinter <- Leo[which(Leo$cat=="Wintering"),]
leoWinter <- leoWinter[c(which(year(leoWinter$timestamp)==2008 &month(leoWinter$timestamp)%in%c(11,12)),
             which(year(leoWinter$timestamp)==2009 &month(leoWinter$timestamp)%in%c(1,2,3))),]
leoWinterNSD <- (spDistsN1(coordinates(leoWinter), coordinates(leoWinter)[1,],longlat=T))^2
plot(leoWinter$timestamp,leoWinterNSD, type="l",
     xlab="Time", ylab="Net square distance (Km²)",main="Winter 2008/2009")
layout(matrix(c(1), ncol=1, byrow=T))
# Note to plot:
# - "All data": shows nicely migration and winter/summer range
# - "Breeding": does not move all to much
# - "Winter08/09": shows that data is not classified correctly, data should be cut off at both ends



#################################
### FIRST PASSAGE TIME (FPT) ####
#################################
library(adehabitatLT)
## project and center the data set to minimize the effects of spherical distortion
Leoprj <- spTransform(Leo, center=T) # center=T: the center of the coordinate system is the center of the track. Units are in meters
# calculate FPT, for different radii, in this case from 1000m to 10^6m in 150 steps
fptLeo <- fpt(as(Leoprj, "ltraj"), radii=10^seq(3, 6, length.out=150), units="days")
# calculate mean nb of days to leave each radii 
meanFPT <- colMeans(fptLeo[[1]], na.rm=T)
radiiFPT <- attributes(fptLeo)$radii
plot(meanFPT~radiiFPT,
     type="l", lwd=2, xlab="Radii in meters",
     ylab="First passage time in days", log="xy")
# Note to plot:
# - with increasing radii size, on average it takes the animal longer to leave the circle
# - we are interested in the changes of slope


### variance of the log(FPT) ######
vars <- varlogfpt(fptLeo, graph=F)
plot(as.numeric(vars)~radiiFPT,
     type="l", lwd=1, lty=2, 
     log="x", ylab="Variance of log first passage time", 
     xlab="Radius in meters")
# Note to plot:
# - minima and maxima can indicate change in the movement process


### fitting LM to min/max peaks of variance of log(fpt)
plot(log10(meanFPT)~log10(radiiFPT),
     type="l", lwd=2, xlab="Log radii in meters",
     ylab="Log first passage time in days")
# fit a model to the largest valley, and largest peak of variance
lm1 <- lm(log10(meanFPT[1:which.min(vars[1:which.max(vars)])])~
            log10(radiiFPT[1:which.min(vars[1:which.max(vars)])]))
lm2 <- lm(log10(meanFPT[which.min(vars[1:which.max(vars)]):which.max(vars)])~
            log10(radiiFPT[which.min(vars[1:which.max(vars)]):which.max(vars)]))
abline(lm1, lty=2)
abline(lm2, lty=3)
text(4, 0.1, paste(signif(summary(lm1)$coefficients[2,1], 2), 
                   "±", 
                   signif(summary(lm1)$coefficients[2,2], 2)), pos=4, cex=0.75)
text(4, 1, paste(signif(summary(lm2)$coefficients[2,1], 2), 
                 "±", 
                 signif(summary(lm2)$coefficients[2,2], 2)), pos=4, cex=0.75)
# Note to plot:
# - flat slope (below 2): directional movement
# - steep slope (around 2): brownian movement


### breaks in the trend of the variance of log(fpt) ####
plot(as.numeric(vars)~radiiFPT,
     type="l", lwd=1, lty=2, 
     ylab="Variance of log first passage time", 
     xlab="Radius in meters", log="x")
breaks <- which(diff(floor(diff(as.numeric(vars))))==-1)+1
abline(v=radiiFPT[breaks])
# Note to plot:
# - besides the local minimum and the global maximum, there is a post maximum hump in the variance


### fitting LM to all changes in slope of variance of log(fpt) ####
plot(log10(meanFPT)~log10(radiiFPT),
     type="n", lwd=4, xlab="Log radii in meters",
     ylab="Log first passage time in days")

lm1 <- lm(log10(meanFPT[1:breaks[1]])~log10(radiiFPT[1:breaks[1]]))
lm2 <- lm(log10(meanFPT[breaks[1]:breaks[2]])~log10(radiiFPT[breaks[1]:breaks[2]]))
lm3 <- lm(log10(meanFPT[breaks[2]:breaks[3]])~log10(radiiFPT[breaks[2]:breaks[3]]))
lm4 <- lm(log10(meanFPT[breaks[3]:breaks[4]])~log10(radiiFPT[breaks[3]:breaks[4]]))
lm5 <- lm(log10(meanFPT[breaks[4]:length(as.numeric(vars))])~log10(radiiFPT[breaks[4]:length(as.numeric(vars))]))

abline(lm1, lty=2, lwd=1 + summary(lm1)$coefficient[2,1], col=alpha("black", 0.8))
abline(lm2, lty=2, lwd=1 + summary(lm2)$coefficient[2,1], col=alpha("black", 0.8))
abline(lm3, lty=2, lwd=1 + summary(lm3)$coefficient[2,1], col=alpha("black", 0.8))
abline(lm4, lty=2, lwd=1 + summary(lm4)$coefficient[2,1], col=alpha("black", 0.8))
abline(lm5, lty=2, lwd=1 + summary(lm5)$coefficient[2,1], col=alpha("black", 0.8))

lines(log10(meanFPT)~log10(radiiFPT),type="l", lwd=4, col=alpha("grey40", 0.8))
legend("bottomright",title="Radii (m)", lty=c(2,2,2,2,2), 
       lwd=signif(c(1+summary(lm1)$coefficient[2,1],
                    1+summary(lm2)$coefficient[2,1],
                    1+summary(lm3)$coefficient[2,1],
                    1+summary(lm4)$coefficient[2,1],
                    1+summary(lm5)$coefficient[2,1]),2),
       c(paste(c(1000, round(radiiFPT[breaks],0))[1:2], collapse=" - "),
         paste(c(1000, round(radiiFPT[breaks],0))[2:3], collapse=" - "),
         paste(c(1000, round(radiiFPT[breaks],0))[3:4], collapse=" - "),
         paste(c(1000, round(radiiFPT[breaks],0))[4:5], collapse=" - "),
         paste(c(round(radiiFPT[breaks],0)[4], 100000), collapse=" - ")),
       bty="n", cex=0.75)
# Note to plot:
# - the different radii represent the different scales at which Leo is operating


### FPT at the 4 different scales ###
par(mfrow=c(2,2))
for(i in 4:1){
  plot(fptLeo[[1]][,breaks[i]]~ Leo$timestamp, type="n",
       xlab="Time", ylab="FPT (days)",
       main=paste("Radius ", round(radiiFPT[breaks[i]],0), "meters"),
       bty="n")
  points(fptLeo[[1]][,breaks[i]]~ Leo$timestamp, pch=16, col=alpha("grey", 0.1))
  lines(fptLeo[[1]][,breaks[i]]~ Leo$timestamp)
}
par(mfrow=c(1,1))
# Note to plot:
# - 4th: ~2Km seems to be Leo's day range size
# - 3rd: at ~20Km, in winter it does not take Leo long to leave 20Km cicle. In summer it takes Leo long time to leave this area, probably this reflects the size of the breeding area of Leo
# - 1st, 2nd: ~70-80Km seems to be the area that Leo occupies in winter



################################################################
## VARIANCE OF dBBMM (dynamic Brownian Bridge Movement Model) ##
################################################################
library(move)
library(lubridate)
Leroy <- move(system.file("extdata","leroy.csv.gz",package="move"))
Leroy <- spTransform(Leroy, center=T)
LeroyVar <- brownian.motion.variance.dyn(Leroy, location.error=25, window.size=71, margin=21)
VarDat <- data.frame(var=getMotionVariance(LeroyVar), hour=hour(LeroyVar$study.local.timestamp))
boxplot(VarDat$var~VarDat$hour, xlab="Hour of the day", ylab="mean Brownian variance", pch="*")
# Note to plot:
# - Leroy moves around during the night, and sleeps during the day



###############################################################
## VARIANCE dBGB (dynamic bi-Gaussian Bridge Movement Model) ##
###############################################################
LeroyBGB <- dynBGBvariance(Leroy, locErr=25, windowSize=31, margin=15)
VarDat <- data.frame(var=getMotionVariance(LeroyBGB), hour=hour(LeroyVar$study.local.timestamp))
VarDat$I_d <- ((VarDat$var.para-VarDat$var.orth)/(VarDat$var.para+VarDat$var.orth))
boxplot(I_d~hour, xlab="Hour of the day", ylab="dBGB variance index", data=VarDat, pch="*")
abline(h=0, lty=2)
# Note to plot:
# - zero: true brownian motion
# - positive values: directional movement, less changes in direction and more in velocity than expected from brownian motion
# - negative values: moving in all directions, there are more changes in direction and less in velocity than expected from brownian motion



##############################
#### SIMULATION OF TRACKS ###
##############################
## simulate Brownian motion
# simulate pure brownian walk
steps <- 1000
duration <- 3600
start.time <- Sys.time()
simmBrown<-   move(x=cumsum(rnorm(steps, 0, 1)),
                   y=cumsum(rnorm(steps, 0, 1)),
                   time=seq(start.time, start.time+duration, length.out=steps),
                   id="SimmBrown")
# simulate biased brownian walk (mean different to 0)
steps <- 1000
duration <- 3600
start.time <- Sys.time()
simmBias<-   move(x=cumsum(rnorm(steps, 0.1, 1)),
                  y=cumsum(rnorm(steps, 0.1, 1)),
                  time=seq(start.time, start.time+duration, length.out=steps),
                  id="SimmBiased")
simsB <- moveStack(list(simmBrown, simmBias))
plot(simsB, type="n")
lines(simsB, col=c("grey", "black"))
legend("topleft", lty=1, col=c("grey", "black"), legend=c("simm. Brownian", "simm. Bias"), bty="n")
# Note to plot:
# - pure brownian walk: has a normal random distribution around a mean of 0, this is a uncorrelated and unbiased random walk, where movement in any direction is equally likely
# - biased walk: mean is different to 0, the walk is directed


### simulate uniform random distribution ###
steps <- 1000
duration <- 3600
start.time <- Sys.time()
simmU <- move(x=cumsum(runif(steps, min=-1.96, max=1.96)),
              y=cumsum(runif(steps, min=-1.96, max=1.96)),
              time=seq(start.time, start.time+duration, length.out=steps),
              id="Simm")
plot(simmU, type="l", xlab="Longitude", ylab="Latitude")
# Note to plot:
# - looks like a brownian random walk, but it isn't. This type of walk is never observed, as normally there is not the same probability to do all types of movement


### comparison of distribution of speeds ###
hist(speed(simmU), col="grey", xlim=c(0,1.1), main=NA, xlab="Speed")
hist(speed(simsB[[1]]), breaks="FD", ylim=c(0,250), add=T, col=alpha("white", 0.5))
legend("topright",fill=c("white","grey"), legend=c("simm. Brownian", "simm. Uniform"), bty="n")
# Note to plot:
# - uniform random distrib.: speeds are skewed to the right, suggesting that the animal is travelling mostly at high speeds => this is highly unnatural
# - normally speeds are left skewed


#### simulate multivariate normal distribution ####
## 1st, create multivariate normal distributed numbers with defined variance/covariance structure, taken from Leo
library(mvtnorm)
sig <- cov(cbind(log(turn$Vm), log(abs(rad(turn$angle)))), use="complete.obs")
mV <- mean(log(turn$Vm), na.rm=T)
mT <- mean(log(abs(rad(turn$angle))), na.rm=T)
rNum <- rmvnorm(10000, mean=c(mV, mT), sigma=sig)
dists <- exp(rNum[,1])
turnVec <- deg(exp(rNum[,2]))%%180
rNeg <- sample(1:length(turnVec), length(turnVec)*0.5)
turnVec[rNeg] <- turnVec[rNeg] * -1
plot(turnVec,dists, pch=16, col=alpha("grey", 0.5))
# Note to plot:
# - function provides realistic looking speeds and turning angles


## 2nd, use the simulated speeds and turning angles to simulate a trajectory 
library(geosphere)
AZ <- (cumsum(c(176, turnVec))%%360)[-1]
start <- matrix(c(0,0), ncol=2)
track <- NULL
for(i in 1:length(AZ)){
  tmp <- destPoint(start, AZ[i], dists[i])
  track <- rbind(track, tmp)
  start <- tmp
}
plot(track, type="l", ylab="Latitude", xlab="Longitude", asp=1)
# Note to plot: 
# - the reason why this track is much less directional than Leo's original, is that although relation between speeds and turning angles is more realistic, the random numbers lack autocorrelation


### simulate correlated random walk ###
library("adehabitatLT")
set.seed(5323)
crw <- simm.crw(1:100, r=.99)

plot(crw)
text(100,30,paste("Correlation coeffcient of azimuth: ", 
                 round(cor(crw[[1]]$abs.angle[-99:-100],
                           crw[[1]]$abs.angle[c(-1,-100)]), 2), 
                 sep=""), pos=2, cex=0.75)
text(100,20,paste("Correlation coeffcient of turning angle: ", 
                 round(cor(crw[[1]]$rel.angle[c(-1, -99, -100)], 
                           crw[[1]]$rel.angle[c(-1, -2 ,-100)]), 2), 
                 sep=""), pos=2, cex=0.75)
# Note to plot:
# - direction is correlated (r=.99), there is a high consistency in direction
# - turning angles are uncorrelated



### randomized track ###
# this function creates random tracks that have the same start and end point as original track, and shuffles the sequence of the original segments, maintaining step lengths and headings in a planar coordinate system but not turning angles #
rand.seg <- function(x, repeats=1){
  coordinateDiffs <- cbind(diff(coordinates(x)[,1]-coordinates(x)[1,1]) , diff(coordinates(x)[,2]-coordinates(x)[1,2]))
  t <- replicate(repeats, apply(rbind(coordinates(x)[1,], coordinateDiffs[sample(1:(n.locs(x)-1)),]), 2, cumsum), simplify=F)
  options(warn=-1)
  MO <- moveStack(lapply(t, function(s) move(x=s[,1], y=s[,2], time=timestamps(x), data=as.data.frame(s), animal="RandTrack")))
  options(warn=0)
  return(MO)
}

library(adehabitatLT)
library(move)
library(scales)
set.seed(43597)
rw <- simm.crw(1:50, h=1, r=0.95, c(10,10))
rw <- as(rw, "Move")
limes <- extent(bbox(rw))*1.3
plot(rw, xlim=c(limes@xmin, limes@xmax), ylim=c(limes@ymin, limes@ymax),
     type="l", lwd=1.6, bty="L",xlab="X-coordinates", ylab="Y-coordinates")

repl <- rand.seg(rw, 99)
lines(repl, col=alpha("black", 0.30), lty=2)
points(10,10, pch=17, cex=1.2)
points(tail(coordinates(rw),1), pch=19)
legend("topleft", pch=c(NA, NA, 17, 19), lty=c(1,2,NA,NA), 
       c("Empirical track", "Randomised track", "Start", "End"), 
       bty="n", cex=0.75)



#################################
## AUTO-CORRELATION STRUCTURE ###
#################################
## correlation between any 2 subsequent speeds (autocorrelation of lag 1)
cor(speed(Leroy)[-length(speed(Leroy))], speed(Leroy)[-1])

## correlation of all speeds separated by one segment (autocorrelation of lag 2)
lag <- 2
cor(speed(Leroy)[seq(1, length(speed(Leroy))-lag, 1)],
    speed(Leroy)[seq(1+lag, length(speed(Leroy)), 1)])
# Note: autocorrelation of lag 2 is lower than of lag 1

## autocorrelation of speeds between lag 1 to lag 100
r <- NULL
p <- NULL
for(lag in 1:100){
  r <- c(r, cor.test(speed(Leroy)[seq(1, length(speed(Leroy))-lag, 1)],
                     speed(Leroy)[seq(1+lag, length(speed(Leroy)), 1)], method="spearman")$estimate)
  p <- c(p, cor.test(speed(Leroy)[seq(1, length(speed(Leroy))-lag, 1)],
                     speed(Leroy)[seq(1+lag, length(speed(Leroy)), 1)], method="spearman")$p.value)
}
plot(r, type="n", xlab="Lag", ylab="Correlation coefficient")
points(r[p<=0.05]~seq(1:100)[p<=0.05], pch=16, col="grey40")
points(r, type="b")
points(40,0.5, pch=16, col="grey40")
points(40,0.5)
text(41,0.49, expression(p<=0.05), pos=4)
# Note to plot:
# - p-value: if correlation coefficient is signif. different from 0
# - significant autocorrelation up to ~lag 10. Corresponds to ~6h
# - oscilation correspond to activity, Leroy is active for a while, and than inactive for a while

## different plot, same information
acf(speed(Leroy))


## autocorrelation of azimuth and turning angle
par(mfcol=c(1,2))
# angular measurements need to be adjusted to their circular nature first
az0 <- as.circular(angle(Leroy), type="direction", 
                   units="degrees", template="geographic", 
                   zero="0", rotation="clock", modulo="asis")
r <- NULL
p <- NULL
for(lag in 1:100){
  r <- c(r, cor.circular(az0[seq(1, length(az0)-lag, 1)],
                         az0[seq(1+lag, length(az0), 1)]))
  p <- c(p, cor.circular(az0[seq(1, length(az0)-lag, 1)],
                         az0[seq(1+lag, length(az0), 1)], test=T)$p.value)
}
plot(r, type="n", xlab="Lag", ylab="Correlation coefficient", ylim=c(-0.1, 1), main="Azimuth")
points(r[p<=0.05]~seq(1:100)[p<=0.05], pch=16, col="grey40")
points(r, type="b")
points(40,0.5, pch=16, col="grey40")
points(40,0.5)
text(41,0.49, expression(p<=0.05), pos=4)

Leroy <- spTransform(Leroy, CRS("+proj=longlat"))
turn0 <- as.circular(turnAngleGc(Leroy), type="angle", 
                     units="degrees", template="geographic", 
                     zero="0", rotation="clock", modulo="asis")
r <- NULL
p <- NULL
for(lag in 1:100){
  r <- c(r, cor.circular(turn0[seq(1, length(turn0)-lag, 1)],
                         turn0[seq(1+lag, length(turn0), 1)]))
  p <- c(p, cor.circular(turn0[seq(1, length(turn0)-lag, 1)],
                         turn0[seq(1+lag, length(turn0), 1)], test=T)$p.value)
}
plot(r, type="n", xlab="Lag", ylab="Correlation coefficient", 
     ylim=c(-0.1, 1), main="Turning angles")
points(r[p<=0.05]~seq(1:100)[p<=0.05], pch=16, col="grey40")
points(r, type="b")
points(40,0.5, pch=16, col="grey40")
points(40,0.5)
text(41,0.49, expression(p<=0.05), pos=4)
# Note to plot:
# - there is little autocorrelation, Leroy moves randomly in all directions


#### Interpolation of positions 
library(scales)
# searching for gaps with missed fixes
maxBreaks <- floor(timeLag(Leroy, units="mins")[timeLag(Leroy, units="mins")>17] / min(timeLag(Leroy, units="mins")[timeLag(Leroy, units="mins")<17]))
minBreaks <- ceiling(timeLag(Leroy, units="mins")[timeLag(Leroy, units="mins")>17] / max(timeLag(Leroy, units="mins")[timeLag(Leroy, units="mins")<17]))
minBreaks[minBreaks==1] <- 2
breaks <- mapply("[", mapply(seq, minBreaks, maxBreaks, 1), lapply(lapply(mapply(seq, minBreaks, maxBreaks, 1), length), sample, size=1))
lags <- (timeLag(Leroy, units="secs")[timeLag(Leroy, units="mins")>17]) / breaks
mts <- lapply(mapply(rep, lags, (breaks-1)), cumsum)
TS <- as.POSIXct(unlist(mapply("+", timestamps(Leroy)[which(timeLag(Leroy, units="mins")>17)], mts)), origin=origin, tz="UTC")
# interpolate fixes at the new timestamps
LeroyInt <- interpolateTime(Leroy, TS, spaceMethod="greatcircle")
LeroyII <- data.frame(x=c(coordinates(Leroy)[,1], coordinates(LeroyInt)[,1]),
                      y=c(coordinates(Leroy)[,2], coordinates(LeroyInt)[,2]),
                      timestamp=c(timestamps(Leroy), timestamps(LeroyInt)))
LeroyII <- LeroyII[order(LeroyII$timestamp),]
LeroyII <- move(x=LeroyII$x, y=LeroyII$y, time=LeroyII$timestamp, id="LeroyInt")
par(mfcol=c(1,2))
plot(Leroy, type="l", col="grey40")
points(Leroy, pch=16, col=alpha("grey", 0.3))
points(LeroyInt, pch=16, col=alpha("black", 0.5), cex=0.5)
legend("topleft",pch=16, legend=c("true","interpolated"),col=c("grey","black"), cex=0.5)
acf(speed(LeroyII), lag.max=100)
# Note to plots:
# - gaps are biologically meaningful: Leroy is denning. By interpolating we are forcing Leroy to move
# - pattern in autocorrelation plot is result from the introduction of these positions

### Distance and time lag between consecutive locations
par(mfrow=c(1,1))
plot(distance(Leroy)[timeLag(Leroy, units="mins")>17],
     timeLag(Leroy)[timeLag(Leroy, units="mins")>17],
     xlab="Distance in meters",
     ylab="Time lag in minutes", log="xy")
# Note to plot:
# - no consistency in speed
# - short distances with large time lags
# - long distances with short time lags
# - if there was a linear relationship, interpolation would make some sense

