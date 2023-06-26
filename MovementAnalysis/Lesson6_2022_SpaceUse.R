#########################################
###           AniMove 2022            ###    
### Script by Kami Safi & Anne Scharf ###
#########################################

library(move)

setwd("/home/ascharf/Documents/Animove22/MovementAnalysis/data/")
bats <- move("Parti-colored bat Safi Switzerland.csv")

###########
## MCP ###
###########
library(adehabitatHR)
X330 <- bats[["X330"]]
X330$id <- "X330"
mcpX330<-mcp(as(X330[,'id'], 'SpatialPointsDataFrame'))
plot(X330, type="n", bty="na", xlab="Longitude", ylab="Latitude")
plot(mcpX330, col="grey90", lty=2, lwd=1.25, add=TRUE)
points(X330, pch=16)
points(X330, pch=1, col="white")
legend("topright", as.character("95% MCP"), fill="grey90", bty="n")
mcpX330
# Note: area value seems strange. That is because our used locations are in the geographic coordinates system (long/lat). adehabitatHR calculated the area according to the units of the projection, in this case decimal degrees


# therefore we have to project our data into a equidistant projection
library("rgeos")
bats$id <- trackId(bats) 
mcpData<-mcp(as(bats[,'id'],'SpatialPointsDataFrame'))
#first option: reproject locations, than calculate mcp
bats.proj <- spTransform(bats, CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +x_0=600000 +y_0=200000 +ellps=bessel +units=m +no_defs"))
mcpData.proj <- mcp(as(bats.proj[, 'id'],'SpatialPointsDataFrame'))
#second option: calculate mcp, than reproject mcp
projection(mcpData) <- CRS("+proj=longlat +datum=WGS84")
mcpData <- spTransform(mcpData, CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +x_0=600000 +y_0=200000 +ellps=bessel +units=m +no_defs"))
plot(bats.proj[["X21"]], bty="na", xlab="Longitude", ylab="Latitude")
plot(mcpData.proj[mcpData.proj$id=="X21",], add=TRUE)
plot(mcpData[mcpData$id=="X21",], add=TRUE, lty=2)
legend("bottomleft", c("First reproject then mcp", "First mcp then reproject"), lty=c(1,2), bty="n")
legend("topleft", sprintf("Area = %.2f", c(gArea(mcpData.proj, byid=TRUE)["X21"],gArea(mcpData, byid=TRUE)["X21"])/1000^2), lty=c(1,2), bty="n")
# Note to plot:
# - the 2 options result in different area calculations
# - adehabitatHR uses distance between locations to do the calculation
# - therefore always project locations, and than calculate MCP


## Size of MCP changes with sampling effort or sampling size
hrBootstrap(bats[['X21']], rep=500, levelMax=95)
legend("bottomright", legend=c("real MCP size","100% percentil","75% percentil","50% percentil","25% percentil","0% percentil"), lty=c(4,2,3,1,3,2), col=c("black","cyan","red","black","red","cyan"))
# Note to plot:
# - if sampling is large enough a saturation between sample size and area is reached



############
## Kernel ##
############
# creating a very simple density plot
library(raster)
template <- raster(extent(bats.proj[[1]]))
res(template)<-500
count <- rasterize(split(bats.proj)[[1]], template,field=1,  fun="count")
plot(count, col=grey(10:0/12))
plot(mcpData.proj[1,], add=TRUE)
points(bats.proj[[1]], pch=16, cex=0.5)
# Note: the result is highly dependent on the chosen cell size. In this case 500x500m


# kernel implementation by "adehabitatHR" library
library(adehabitatHR)
library(scales)
X21 <- bats.proj[['X21']]
kern1 <- kernelUD(as(X21, "SpatialPoints"), h=500)
kern2 <- kernelUD(as(X21, "SpatialPoints"))
kern3 <- kernelUD(as(X21, "SpatialPoints"), h=2000)
kern4 <- kernelUD(as(X21, "SpatialPoints"), h="LSCV")
par(mfrow=c(2,2))
par(mar=c(1,0.5,3,0.5))
kern <- c("kern1", "kern2", "kern3", "kern4")
hName <- c("h=500",
           "h='ad-hoc'",
           "h=2000",
           "h=LSCV")
for(i in 1:4){
  plot(getverticeshr(get(kern[i])))
  points(X21, pch=16, cex=0.75, col=alpha("black", 0.2))
  points(X21, cex=0.75)
  title(hName[i])
}
# Note to plot:
# - h: degree of smoothness or how tightly the data should be hugged by the distribution function
# - h="LSCV": h calculated from the data via least square cross validation
# - h="ad-hoc": h calculated from the data via sample size and spatial spread


# kernel implementation by "ks" library
library(ks) 
library(scales)
pos <- coordinates(X21)
H.bcv <- Hbcv(x=pos)
H.pi <- Hpi(x=pos) 
H.lscv <- Hlscv(x=pos) 
H.scv <- Hscv(x=pos) 

par(mfrow=c(2,2))
par(mar=c(1,0.5,3,0.5))
H <- c("H.bcv", "H.pi", "H.lscv", "H.scv")
hT <- c("Biased cross-validation (BCV)",
        "Plug-in",
        "Least-squares cross-validation",
        "Smoothed cross-validation")
for(i in 1:4){
  fhat <- kde(x=pos, H=get(H[i]), compute.cont=TRUE) 
  plot(fhat, cont=c(75, 50, 5), bty="n", 
       xaxt="n", yaxt="n", 
       xlab=NA, ylab=NA, asp=1,display="filled.contour")
  points(X21, pch=16, cex=0.75, col=alpha("black", 0.2))
  title(hT[i])
}



###########
## LoCoH ##
###########
# check vignettes
library(move)
library(maptools)
library(adehabitatHR)
data(leroy)
# data need to be transformed into a equidistant projection because method relays on distance calculations
leroy <- spTransform(leroy, center=TRUE)

# png(filename="locoh_maps.png", width = 21, height = 15, units = "cm", res=300)
par(list(mfrow=c(2,2), mar=c(2,2,2,2)))
leroy.mcp <- mcp(as(leroy, "SpatialPoints"), percent=95)
plot(leroy.mcp, col=grey(0.9), lty=2, lwd=2)
points(leroy, col="#00000060", pch=16, cex=0.5)
lines(leroy, col="#00000030")
title("Minimum convex polygon")
# include "k" number of closest neighbour locations
kLoc <- LoCoH.k(as(leroy, "SpatialPoints"), k=75)
plot(kLoc, col=grey((0:length(kLoc)/length(kLoc))*0.7), border=NA)
title("k-NNCH LoCoH")
# include location within a radius "r"
rLoc <- LoCoH.r(as(leroy, "SpatialPoints"), r=800)
plot(rLoc, col=grey((0:length(rLoc)/length(rLoc))*0.7), border=NA)
title("r-NNCH LoCoH")
# sum of distances of included neighbour locations to root location is "a"
aLoc <- LoCoH.a(as(leroy, "SpatialPoints"), a=9000)
plot(aLoc, col=grey((0:length(aLoc)/length(aLoc))*0.7), border=NA)
title("a-NNCH LoCoH")
# dev.off()

library(imager)
dev.off()
plot(load.image("locoh_maps.png"),axes=F)

## area changes depending on the choice of k, r, or a
dev.off()# to reset the plotting environment
# png(filename="locoh_Area.png", width = 21, height = 15, units = "cm", res=300)
par(mfrow=c(1,3))
kLocArea <- LoCoH.k.area(as(leroy, "SpatialPoints"), 
                         krange=floor(seq(75, 500, length=10)), 
                         percent=90)
title("k-NNCH LoCoH")
rLocArea <- LoCoH.r.area(as(leroy, "SpatialPoints"), 
                         rrange=seq(500, 1600, 100), 
                         percent=90)
title("r-NNCH LoCoH")
aLocArea <- LoCoH.a.area(as(leroy, "SpatialPoints"), 
                         arange=seq(5000, 13000, 1000), 
                         percent=90)
title("a-NNCH LoCoH")
# dev.off()

plot(load.image("locoh_Area.png"),axes=F)

##############
## t- LoCoH ##
##############
# check website and vignettes
# for installation go to: http://tlocoh.r-forge.r-project.org/#installation
library(tlocoh)
leroy.lxy <- move.lxy(leroy) # tlocoh has its own object, class lxy

## 1. without including time
## calculate a series of hullsets with different number of "k" nearest neighbours. by setting s=0, time is ignored
leroy.lxy <- lxy.nn.add(leroy.lxy, s=0, k=seq(5, 105, 10))
leroy.lhs <- lxy.lhs(leroy.lxy, k=seq(5, 105, 10), s=0)
leroy.lhs <- lhs.iso.add(leroy.lhs)


## plotting the results from above to see how the value of "k" affects area
# png(filename="tLoCoH_noTime.png", width = 21, height = 15, units = "cm", res=300)
par(mfrow=c(1,2))
par(list(mar=c(5, 4, 4, 2) + 0.1), bty="n")
iso.info.all <- do.call(rbind, lapply(leroy.lhs, function(myhs) do.call(rbind, lapply(myhs$isos, function(myiso) data.frame(id = myhs[["id"]], mode = myhs[["mode"]], s = myhs[["s"]], param.val = myhs[[myhs[["mode"]]]], sort.metric = myiso[["sort.metric"]], myiso[["polys"]]@data[c("iso.level","area")])))))
iso.info.all$area <- iso.info.all$area/(1000*1000)
plot(area~param.val, type="n", data=iso.info.all, xlab="Number of neighbours (k)", ylab=expression(paste("Area in ", km^2, sep="")), ylim=c(-0.1, 16.1))
for(i in 1:length(unique(iso.info.all$iso.level))){
  tmp <- iso.info.all[iso.info.all$iso.level==unique(iso.info.all$iso.level)[i],]
  lines(area~param.val, type="l", data=tmp, lty=i)
}
legend("topleft", as.character(unique((iso.info.all$iso.level))), lty=1:length(unique(iso.info.all$iso.level)), cex=0.6, bty="n")
par(mar=c(1,1,1,1))
plot(leroy.lhs[[8]]$isos[[1]]$polys, col=grey((0.7*seq(0.1,0.99,length=5))), border=NA)
points(leroy, col="#00000060", pch=16, cex=0.5)
lines(leroy, col="#00000030")
text(0,3000,"k-NNCH LoCoH for k=75", pos=4, cex=0.75)
legend(-1636, -2018, as.character(unique((iso.info.all$iso.level))), fill=grey((0.7*seq(0.1,0.99,length=5))), cex=0.6, bty="n")
# dev.off()
plot(load.image("tLoCoH_noTime.png"),axes=F)

## the calculation for hullset with "a" cumulative distance to nearest neighbours, is done the same as above, but instead of argument "k", we use argument "a". 

## 2. including time
## plot to decide which s value to choose, based on a temporal scale
# png(filename="tLoCoH_selectTime.png", width = 21, height = 15, units = "cm", res=300)
leroy.lxy <- move.lxy(leroy)
leroy.lxy <- lxy.ptsh.add(leroy.lxy)
# dev.off()
plot(load.image("tLoCoH_selectTime.png"),axes=F)
# Note: select a "s" so that 40-60% of the hulls are time selected. In this case ~ 0.03
# ptsh:proportion of time-selected hulls
# "s" depends on the sampling schedule and the map units


## plot to decide which s value to choose, based on time scaled distance
# png(filename="tLoCoH_selectTime2.png", width = 21, height = 15, units = "cm", res=300)
lxy.plot.sfinder(leroy.lxy, delta.t=3600*c(6,12,24,36,48,54,60))
# dev.off()
plot(load.image("tLoCoH_selectTime2.png"),axes=F)
# Note: s=0.03 kind of fits with the daily behaviour


## again calculate a series of hullsets with different number of "k" nearest neighbours, but this time with the estimated "s" value
leroy.lxy <- lxy.nn.add(leroy.lxy, s=0.03, k=seq(10, 100, 10))
leroy.lhs <- lxy.lhs(leroy.lxy, s=0.03, k=seq(10, 100, 10))
leroy.lhs <- lhs.iso.add(leroy.lhs, k=seq(10, 100, 10))

## plot the results from above
lhs.plot.isoear(leroy.lhs)
# Note: the amount of edge (holes) should not be to large. Probably want to choose 1st minima
plot(leroy.lhs, iso=TRUE, record=TRUE)
# Note (both plots): isopleth 10% contains 10% of locations; isopleth 100% contains 100% of the locations


## different plot but containing same info as above
# png(filename="tLoCoH_Time.png", width = 21, height = 15, units = "cm", res=300)
par(mfrow=c(1,2))
par(list(mar=c(5, 6, 4, 2) + 0.1), bty="n")
iso.info.all <- do.call(rbind, lapply(leroy.lhs, function(myhs) do.call(rbind, 
        lapply(myhs$isos, function(myiso) data.frame(id = myhs[["id"]], 
            mode = myhs[["mode"]], s = myhs[["s"]], param.val = myhs[[myhs[["mode"]]]], 
            sort.metric = myiso[["sort.metric"]], myiso[["polys"]]@data[c("iso.level", "area", "edge.len")])))))
plot(I(edge.len/area)~param.val, type="n", data=iso.info.all, 
     xlab="Number of neighbours (k)", 
     ylab=expression(edge%/%area))
for(i in 1:length(unique(iso.info.all$iso.level))){
  tmp <- iso.info.all[iso.info.all$iso.level==unique(iso.info.all$iso.level)[i],]
  lines(I(edge.len/area)~param.val, type="l", data=tmp, lty=i)
}
legend("topright", as.character(unique((iso.info.all$iso.level))), lty=1:length(unique(iso.info.all$iso.level)), cex=0.6, bty="n")
par(mar=c(1,1,1,1))
plot(leroy.lhs[[5]]$isos[[1]]$polys, col=grey((0.7*seq(0.1,0.99,length=5))), border=NA)
points(leroy, col="#00000060", pch=16, cex=0.5)
lines(leroy, col="#00000030")
text(0,3000,"k-NNCH LoCoH for k=50", pos=4, cex=0.75)
legend(-1636, -2018, as.character(unique((iso.info.all$iso.level))), fill=grey((0.7*seq(0.1,0.99,length=5))), cex=0.6, bty="n")
par(mfrow=c(1,1))
# dev.off()
plot(load.image("tLoCoH_Time.png"),axes=F)

################
## dBBMM & UD ##
################
library(move)
leroy <- move(system.file("extdata","leroy.csv.gz",package="move"))
leroy <- spTransform(leroy, center=TRUE)
# check timeLag of data. If there are timelags shorter than intended, use the argument "timestep" in the dBBMM function (see ?brownian.bridge.dyn), to make sure calculation does not take forever
summary(timeLag(leroy,"mins")) 
BB.leroy <- brownian.bridge.dyn(leroy, ext=.45, dimSize=150, location.error=20, margin=11, window.size=31)
plot(BB.leroy)
# Note to plot: 
# - all pixels sum up to 1 (total time tracked)
# - the values correspond to the proportion of the time tracked that the animal spend in that area
# - in this case, a pixel with value 0.035 means that leroy spend 3.5% from the total tracking time in that pixel 
# - these results are very useful to find out where the animal spend how much time, and e.g. relate it to environmental variables

## extract the utilization distribution (UD)
udleroy <- getVolumeUD(BB.leroy)
plot(udleroy, col=terrain.colors(100))

## from the ud object, also the contours can be extracted
plot(leroy, col="#00000060", pch=16, cex=0.5, bty="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA)
lines(leroy, col="#00000030")
contour(udleroy, levels=c(0.5, 0.95), add=TRUE, lwd=c(2, 1), lty=c(2,1))
title("Dynamic brownian bridge")

## plotting the UD95 on google map
library(ggmap)
library(ggplot2)
cl95 <- raster2contour(BB.leroy, levels=0.95)
cl95LL <- spTransform(cl95, CRS("+proj=longlat"))
cl95df <- data.frame(do.call("rbind", coordinates(cl95LL)[[1]]))
cl95df$polyN <- rep(1:length(coordinates(cl95LL)[[1]]), lapply(coordinates(cl95LL)[[1]], nrow))
leroyDF <- as.data.frame(spTransform(leroy,"+proj=longlat"))
m <- get_map(sp::bbox(extent(cl95LL)*1.5), zoom=13, source="google", maptype="hybrid")
ggmap(m)+geom_path(data=cl95df, aes(x=X1,y=X2,group=polyN),color="red")+
  geom_path(data=leroyDF, aes(x=location.long, y=location.lat),alpha=0.2)+
  geom_point(data=leroyDF, aes(x=location.long, y=location.lat),alpha=0.3, shape=20)+
  labs(x="",y="")+
  theme(axis.text=element_blank(),axis.ticks=element_blank())+ 
  theme(legend.position="none")



## checking for effect of margin and window size on the results
par(mfrow=c(3,3), mar=c(1,2,2,1))
margins <- c(15, 9, 3)
windows <- c(101, 67, 33)
runs <- expand.grid(margins, windows)
for(i in 1:nrow(runs)){
  BB.leroy <- brownian.bridge.dyn(leroy, dimSize=150, location.error=20, margin=runs[i,1],
                                  window.size=runs[i,2], time.step=2, ext=2)
  udleroy <- getVolumeUD(BB.leroy)
  udleroy99 <- udleroy
  udleroy99[udleroy99>0.99] <- NA
  udleroy99 <- trim(udleroy99,padding=5)
  contour(udleroy99, levels=c(0.5, 0.95), bty="n", xaxt="n", 
          yaxt="n", xlab=NA, ylab=NA, asp=1)
  mtext(paste("Margin = ", runs[i,1], sep=""), 2)
  mtext(paste("Window size = ", runs[i,2]), 3)
}
# Note to plot: no need to worry all to much about the window size and the margin, as they do not have a major impact on the results. Default values work well.



### effect of sampling frequency and tracking duration on area
dev.off()
par(mfrow=c(1,1))
par(list(mar=c(5, 4, 4, 2) + 0.1), bty="o")
set.seed(3628492)
steps <- 100000
prop=seq(0.01,1,0.01)

track <- as(simm.crw(date=1:steps, h=1, r=0.8), "Move")
thin <- lapply(prop, function(x) track[round(seq(1, steps, length.out=steps * x)), ])
short <- lapply(prop, function(x) track[1:round(steps * x),])
ThinAreas <- lapply(lapply(lapply(thin, as, "SpatialPoints"), kernelUD), kernel.area, percent=95)
ShortAreas <- lapply(lapply(lapply(short, as, "SpatialPoints"), kernelUD), kernel.area, percent=95)

plot(I(unlist(ThinAreas)/min(unlist(ThinAreas)))~seq(1:100), 
     xlab="Percent of the track", 
     ylab="Relative area", ylim=c(0,1.75), 
     type="l", lwd=2, lty=1)
lines(I(unlist(ShortAreas)/max(unlist(ShortAreas)))~seq(1:100), lty=2, lwd=2)
abline(h=1)
legend("topright", c("Thinned trajectory", "Shortened trajectory"), lty=c(1,2), lwd=c(2,2), bty="n")
# Note to plot: 
# - in this case using the kernelUD, when sampling frequency is lower, i.e. thinned trajectory, than the estimated UD is larger than when the sampling frequency is higher
# - the longer the trajectory, the larger the UD


