# development branch of ctmm (more recent than CRAN)
remotes::install_github("ctmm-initiative/ctmm")
# or
devtools::install_github("ctmm-initiative/ctmm")

# ctmm user group for any questions or help
browseURL("https://groups.google.com/g/ctmm-user")

# ctmm point-and-click app - if you know anyone that doesn't user R
devtools::install_github("ctmm-initiative/ctmmweb")
ctmmweb::app()

#! load the ctmm package
library(ctmm)

# load buffalo data from Movebank CSV (which can be compressed)
Leo <- as.telemetry('Leo-65545.csv.gz')
# you can also import from a move object, data.frame, etc.

#! load buffalo dataset from ctmm
data(buffalo)

# this is a list of buffalo telemetry objects
class(buffalo)

# number of buffalo datasets
length(buffalo)

# names of buffalo
names(buffalo)

# summarize buffalo data
summary(buffalo)

###################
# PLOT TELEMETRY
###################

# plot all buffalo
plot(buffalo)
# but all the same color

# plot buffalo with list-sorted rainbow of colors
COL <- rainbow(length(buffalo))
plot(buffalo,col=COL)

# plot buffalo with spatially-separated rainbow of colors
COL <- color(buffalo,by='individual')
plot(buffalo,col=COL)

# many other built in coloring options for telemetry objects
?color
# you can color by sunlight, moonlight, season, time, ...

####################
# PROJECTIONS
####################

# what projection are the buffalo in
projection(buffalo)
compass()

# You want a projection that is locally flat over your data (to minimize distortion).
# By default, as.telemetry() will choose a two-point equidistant projection, which is
# safer for migratory species, but does not preserve North=up.

# center the projection on the geometric median of the data
projection(buffalo) <- median(buffalo)

projection(buffalo)

# now North=up, which is fine for this dataset
plot(buffalo)
compass()

####################
# OUTLIER DETECTION
####################
# Identify outliers consistent/inconsistent with an error model

#! load turtle data
data(turtle)

# this is a list of turtle (and calibration) data
class(turtle)

# number of datasets within list
length(turtle)

# names of each dataset
names(turtle)

#! select the third dataset (female 231)
DATA <- turtle$F231

# plot turtle
plot(DATA)
# notice the error circles

# help file for outlie function
?outlie
# note the 'by' argument in particular

#! calculate outlier statistics and generate plot
OUT <- outlie(DATA)
# red segments are for speed
# blue points are for proximity
# the speed and proximity estimates are error informed (if your data are calibrated)

# outlier statistics (used for coloring)
plot(OUT)
# speed is measured in meters/second
# distance (proximity) is measured in meters

#! index of fix with highest speed
BAD <- which.max(OUT$speed)

#! remove the outlier from the dataset
DATA <- DATA[-BAD,]
# negative indices remove those indices

# plot data with outlier removed
plot(DATA)

# look for outliers again
OUT <- outlie(DATA)
plot(OUT)

###################
# VARIOGRAM
###################

# names of buffalo
names(buffalo)

#! select buffalo Cilla
DATA <- buffalo$Cilla

# plot telemetry object
plot(DATA)

# color by time
COL <- color(DATA,by='time')
plot(DATA,col=COL)

#! calculate a variogram object (named SVF) from the telemetry object
SVF <- variogram(DATA)
plot(SVF)
# on average how far apart (in distance^2) given a time lag between

# help file for variogram
?variogram
# there are some options in here if you have very irregular data
# res, fast, dt
vignette('variogram')
# Sec. "Irregular Sampling Schedules"

# more accurate CIs, slow for larger datasets
SVF <- variogram(DATA,CI="Gauss")

# frequently you want to zoom in to the beginning of the variogram
# plot with zoom slider
zoom(SVF)
# things to look for
# * the asymptote (if any)
# * how long does it take to asymptote
# * initial curvature or initial linear?


###################
# MODEL SELECTION
###################

# model guesstimate function
?ctmm.guess
# variogram will be calculated automatically (with default arguments)
# this is interactive mode
ctmm.guess(DATA)
# notice how much work I spent automating the units of every plot

# this is noninteractive mode
GUESS <- ctmm.guess(DATA,interactive=FALSE)

# automated model selection
?ctmm.select

# fit a bunch of models, tell me what models are being fit, return all models, and use all but one CPU core
FITS <- ctmm.select(DATA,GUESS,trace=3,verbose=TRUE,cores=-1)
# candidate models: OUF, OUf, OUΩ, IOU, BM, IID, inactive
# I've already run this code for you
# save(FITS,file="cilla.rda")
load("cilla.rda")

# lets look at the results
summary(FITS)

# IID was not attempted because the nested-model hierarchy is OUF -> OU -> IID
# so let's include the IID models
FITS[["IID anisotropic"]] <- ctmm.fit(DATA)
FITS[["IID"]] <- ctmm.fit(DATA,ctmm(isotropic=TRUE))

# now including IID model
summary(FITS)

# lets look at individual models
# IID  anisotropic model
summary(FITS[[5]])

# compare mean and covariance to data
plot(DATA,FITS[[5]])

# compare empirical variogram to that of model
zoom(SVF,FITS[[5]])

# calculate residuals
RES <- residuals(DATA,FITS[[5]])

# scatter plot of residuals
plot(RES)

# calculate correlogram of residuals
ACF <- correlogram(RES,res=10)
# res=10 is for drifting sampling rate

zoom(ACF)

# The selected OUF anisotropic model
summary(FITS[[1]])
# area here is Gaussian area
# speed here is Gaussian RMS speed

plot(DATA,FITS[[1]])

zoom(SVF,FITS[[1]])
# not perfect, but much better

RES2 <- residuals(DATA,FITS[[1]])

ACF2 <- correlogram(RES2,res=10)

zoom(ACF2)

# looking at the other models
zoom(SVF,FITS$OUF)

# why is this model fit deflected down?
zoom(SVF,FITS$`OU anisotropic`)


################
# TEASER
################

# simulate data from the selected model with same times
SIM <- simulate(FITS[[1]],t=DATA$t)

# plot data
plot(SIM)
# what areas does this individual like/dislike?

plot(SIM,FITS[[1]],level=NA)

############
# Inês will next cover:
# * AKDE home-range estimation, overlap, meta-analysis
# * experimental design

##########
# Other topics that we can cover (time permitting):
# * range distributions (KDE,AKDE,MCP,...) versus occurrence distribution (Brownian bridge,...)
# * location error modeling & outlier detection
# * speed estimation, diffusion estimation, ...
# * rsf.fit, RSF-AKDEs, AKDEs with boundaries, ...
# * simulations, conditional simulations, predictions, occurrence distributions
# * population distribution estimation
# * exporting parameters + uncertainties for meta-analytic/random-effect regression
# * periodicities


###########
# RANGE VERSUS OCCURRENCE DISTRIBUTIONS
###########

# include Brownian motion models
FITS[["BM"]] <- ctmm.fit(DATA,ctmm(tau=Inf,isotropic=TRUE))
# this one is not as commonly used, but let's throw it in
FITS[["BM anisotropic"]] <- ctmm.fit(DATA,ctmm(tau=Inf))

# you can't compare stationary (IID,OU,OUF) and conditionally stationary (BM,IOU) models with likelihood
summary(FITS)
# but you can compare within
summary(FITS[c("BM","BM anisotropic")])

# again, the selected model looks okay
zoom(SVF,FITS[[1]])

# the Brownian motion model...
zoom(SVF,FITS$BM)
# zoom in

# range distribution - using the selected model
RD <- akde(DATA,FITS[[1]])

# occurrence distribution - using the selected model
OD <- occurrence(DATA,FITS[[1]])

# conventional (non-dynamic) Brownian bridge
BB <- occurrence(DATA,FITS$BM)

# plot them
EXT <- extent(list(DATA,OD,RD))
plot(RD,col.level=NA,col.grid=NA,ext=EXT)
title("OUF Krige")
# plot OUF occurrence distribution
plot(OD,col.level=NA,ext=EXT)
title("OUF Krige")
# plot BM occurrence distribution (BB)
plot(BB,col.level=NA,ext=EXT)
title("BM Krige (BB)")

# Q: What is the occurrence distribution?
# A: Given a random time *in the sampling period*, where was the animal

# Q: What is the range distribution?
# A: At some time in the future/past *under the same behaviors* where will the animal be
# A: Long-term space use *for continuing behaviors*

# Impact of coarsening the data
SUB <- DATA
# remove every other time
SUB <- SUB[as.logical(1:nrow(SUB)%%2),]
par(mfrow=c(1,2))
RD <- akde(SUB,FITS[[1]])
OD <- occurrence(SUB,FITS[[1]])
plot(RD,col.level=NA,col.grid=NA,ext=EXT)
title("Range distribution")
plot(OD,col.level=NA,ext=EXT)
title("Occurrence distribution")

# how much data when they look similar?
nrow(SUB)

# Impact of truncating the data
SUB <- DATA
# remove the second half of the data
SUB <- SUB[1:round(nrow(SUB)/2),]
par(mfrow=c(1,2))
RD <- akde(SUB,FITS[[1]])
OD <- occurrence(SUB,FITS[[1]])
plot(RD,col.level=NA,col.grid=NA,ext=EXT)
title("Range distribution")
plot(OD,col.level=NA,ext=EXT)
title("Occurrence distribution")

par(mfrow=c(1,1))
# range area = predicted space use, given the same behaviors (biological)
# occurrence area = uncertainty (sampling dependent and limited to the sampling period)
# neither estimate space use during the sampling period!!!

##########################
# Occurrence distributions
##########################

library(ctmm)
data(buffalo)
projection(buffalo) <- median(buffalo)
DATA <- buffalo$Cilla
load("cilla.rda")

plot(DATA)

OD <- occurrence(DATA,FITS[[1]])
plot(OD,col.level=NA)

SIM <- simulate(DATA,FITS[[1]],dt=5 %#% 'min')
plot(SIM)

# sum(RASTER*OD) = E[RASTER]
