#########################################
###           AniMove 2022            ###    
### Script by Kami Safi & Anne Scharf ###
#########################################

library(move)
setwd("/home/ascharf/Documents/Animove22/MovementAnalysis/data")

###################################
## GETTING TRACKING DATA INTO R ###
###################################

###--------------------------------------------####
### 1. Directly downloading data from Movebank ####
###--------------------------------------------####

### store the movebank credentials
cred <- movebankLogin(username="RBook", password="Obstberg1")
## or
cred <- movebankLogin() # specially useful when sharing scripts with colleagues

### browse the database ###
## search for studies using keywords included in the study name
searchMovebankStudies(x="bat", login=cred)
searchMovebankStudies(x="Parti-colored bat", login=cred)

# ## if previous function produced an error, than:
# #### ---- set the curlHandle if necessary ---------------------------#######
# curl <- getCurlHandle()
# options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem",
#                                                  package = "RCurl"),
#                             ssl.verifypeer = FALSE, ssl.verifyhost = FALSE))
# curlSetOpt(.opts = list(proxy = 'proxyserver:port'), curl = curl)
# ##### ---------------------------------------------------------------######


## get the metadata of the study
getMovebankStudy(study="Parti-colored bat Safi Switzerland",login=cred)

## check for reference data of animals, deployments and tags
getMovebankReferenceTable(study="Parti-colored bat Safi Switzerland",login=cred)[1:4,]
## check reference data of animals
getMovebankAnimals(study="Parti-colored bat Safi Switzerland",login=cred)[1:4,]

### Download the location data ##
## get the all data
bats <- getMovebankData(study="Parti-colored bat Safi Switzerland", login=cred)
bats

## get only bat "191"
bat191 <- getMovebankData(study="Parti-colored bat Safi Switzerland", login=cred, 
                          animalName="191")
bat191

## get data for a specific time range e.g. between "2002-06-02 23:06:15"
## and "2002-06-11 22:18:25". Time format: 'yyyyMMddHHmmssSSS'
bats.subset <- getMovebankData(study="Parti-colored bat Safi Switzerland",
                               login=cred,
                               timestamp_start="20020602230615000",
                               timestamp_end="20020611221825000")
bats.subset


###----------------------------------------------------###
### 2. Reading in a .csv file downloaded from Movebank ###
###----------------------------------------------------###
bats <- move("Parti-colored bat Safi Switzerland.csv")
bats

## ---- also read EnvData .zip files or tar-compressed csv exports -----------------------------
batsTemp <- move("Parti-colored bat Safi Switzerland-5752797914261819198.zip")
# this data set was annotated with temperature data with the EnvData tool from Movebank

## NOTE: the 'move' function also reads in objects classes from other packages 
# (ltraj, telemetry, track_xyt, track, binClstPath, data.frame). See ?move for more details.


###---------------------------------------------###
### 3. Creating a move object from any data set ###
###---------------------------------------------###
# read the data and store in a data frame
file <- read.csv("Parti-colored bat Safi Switzerland.csv", as.is=T)
str(file)

# first make sure the date/time is in POSIXct format
file$timestamp <- as.POSIXct(file$timestamp,format="%Y-%m-%d %H:%M:%S",tz="UTC")

# also ensure that timestamps and individuals are ordered
file <- file[order(file$individual.local.identifier, file$timestamp),]

# convert a data frame into a move object
Bats <- move(x=file$location.long,y=file$location.lat,
             time=file$timestamp, #already in POSIXct format
             data=file, proj=crs("+proj=longlat +ellps=WGS84"),
             animal=file$individual.local.identifier, sensor="gps")
Bats


###---------------------------------------------###
### download data from Movebank Data Repository ###
###---------------------------------------------###
reposDuck <- getDataRepositoryData("doi:10.5441/001/1.2k536j54")
reposDuck


###---------------------------------------------------###
### example of how to deal with duplicated timestamps ###
###---------------------------------------------------###

### buffalo data 
buffalo <- move("Kruger African Buffalo, GPS tracking, South Africa.csv.gz")

## one solution is to use the argument removeDuplicatedTimestamps=T. But use with care! as duplicates are removed randomly. Additional information could be lost.
buffaloNoDupl <- move("Kruger African Buffalo, GPS tracking, South Africa.csv.gz", 
                      removeDuplicatedTimestamps=T)

## example to remove duplicated timestamps in a controlled way:
## create a data frame
buffalo.df <- read.csv('Kruger African Buffalo, GPS tracking, South Africa.csv.gz', as.is=TRUE)

## get a quick overview
head(buffalo.df, n=2)

## first make sure the date/time is in the correct format
buffalo.df$timestamp <- as.POSIXct(buffalo.df$timestamp, format="%F %T ", tz="UTC")

## also ensure that timestamps and individuals are ordered
buffalo.df <- buffalo.df[order(buffalo.df$individual.local.identifier, buffalo.df$timestamp),]

## get the duplicated timestamps
dup <- getDuplicatedTimestamps(buffalo.df)

## get an overview of the amount of duplicated timestamps
table(unlist(lapply(dup,function(x)length(x)))) 

## inspect the first duplicate
dup[1]
buffalo.df[dup[[1]],] # exact duplicate

## drop exact duplicates taking into account all columns besides "event.id"
buffalo.clean <- buffalo.df[!duplicated(buffalo.df[,!names(buffalo.df) %in% "event.id"]),]

## check again for duplicated timestamps
dup <- getDuplicatedTimestamps(buffalo.clean)
buffalo.clean[dup[[1]],] # same time, different locations

## we will keep the position that results in the shortest distance between the previous and the next location
## A while loop will ensure that the loop continues until each duplicate is removed
## ==> loop starts here
while(length(dup <- getDuplicatedTimestamps(buffalo.clean))>0){
  allrowsTOremove <- lapply(1:length(dup), function(x){
    # row numbers of duplicates
    rown <- dup[[x]]
      # create a row number ID to find the duplicated timestamps in the subset per individual
      buffalo.clean$rowNumber <- 1:nrow(buffalo.clean)
      # subset for the individual, as distances should be measured only within the individual
      ind <- unlist(strsplit(names(dup[x]),split="|", fixed=T))[1]
      subset <- buffalo.clean[buffalo.clean$individual.local.identifier==ind,]
      # if the duplicated positions are in the middle of the table
      if(subset$rowNumber[1]<rown[1] & subset$rowNumber[nrow(subset)]>max(rown)){
        # calculate total distance between fix before/after and the first alternate location
        dist1 <- sum(distHaversine(subset[subset$rowNumber%in%c((rown[1]-1),(max(rown)+1)),c("location.long", "location.lat")],
                                   subset[subset$rowNumber==rown[1],c("location.long", "location.lat")]))
        # calculate total distance between fix before/after and the second alternate location
        dist2 <- sum(distHaversine(subset[subset$rowNumber%in%c((rown[1]-1),(max(rown)+1)),c("location.long", "location.lat")],
                                   subset[subset$rowNumber==rown[2],c("location.long", "location.lat")]))
        # omit the alternate location that produces the longer route
        if(dist1<dist2){rowsTOremove <- rown[2]}else{rowsTOremove <- rown[1]}
      }
      
      # in case the duplicated timestamps are the first location of the animal, calculate distance between duplicate and following location
      if(subset$rowNumber[1]==rown[1]){
        dist1 <- sum(distHaversine(subset[subset$rowNumber==(max(rown)+1),c("location.long", "location.lat")],
                                   subset[subset$rowNumber==rown[1],c("location.long", "location.lat")]))
        dist2 <- sum(distHaversine(subset[subset$rowNumber==(max(rown)+1),c("location.long", "location.lat")],
                                   subset[subset$rowNumber==rown[2],c("location.long", "location.lat")]))
        # and omit the alternate location that produces the longer route
        if(dist1<dist2){rowsTOremove <- rown[2]}else{rowsTOremove <- rown[1]}
      }
      
      # in case the duplicated timestamps are the last positions, calculate distance between duplicate and previous location
      if(subset$rowNumber[nrow(subset)]==max(rown)){
        dist1 <- sum(distHaversine(subset[subset$rowNumber==(rown[1]-1),c("location.long", "location.lat")],
                                   subset[subset$rowNumber==rown[1],c("location.long", "location.lat")]))
        dist2 <- sum(distHaversine(subset[subset$rowNumber==(rown[1]-1),c("location.long", "location.lat")],
                                   subset[subset$rowNumber==rown[2],c("location.long", "location.lat")]))
        # and omit the alternate location that produces the longer route
        if(dist1<dist2){rowsTOremove <- rown[2]}else{rowsTOremove <- rown[1]}
      }
    return(rowsTOremove)
  })
  buffalo.clean <- buffalo.clean[-unique(sort(unlist(allrowsTOremove))),]
  buffalo.clean$rowNumber <- NULL
}
## ==> and ends here

# define the data.frame as a move object after cleaning
buffalo <- move(x=buffalo.clean$location.long,
                y=buffalo.clean$location.lat,
                time=buffalo.clean$timestamp,
                data=buffalo.clean,
                proj=crs("+proj=longlat +datum=WGS84"),
                animal=buffalo.clean$individual.local.identifier,
                sensor=buffalo.clean$sensor.type)



#################
## PROJECTION ###
#################
## check projection 
projection(Bats)
head(coordinates(Bats))

## re-project 
BatsProj <- spTransform(Bats, CRSobj="+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
projection(BatsProj)
head(coordinates(BatsProj))


##############################
### MOVE/MOVESTACK OBJECTS ###
##############################

### Move object (1 indiv) ####
bat191
str(bat191)
## access the information 
timestamps(bat191)[1:5]
coordinates(bat191)[1:5,]
namesIndiv(bat191)
n.locs(bat191)
idData(bat191)
projection(bat191)
extent(bat191)


### Movestack object (multiple indiv) #####
buffalo
str(buffalo)
## access the information 
trackId(buffalo)[1:5]
n.indiv(buffalo)
namesIndiv(buffalo)
n.locs(buffalo)
idData(buffalo)

## get a specific individual (subset is a Move)
Queen <- buffalo[['Queen']]
class(Queen)
## or several (subset is a MoveStack)
CillaGabs <- buffalo[[c("Cilla",'Gabs')]]
class(CillaGabs)

## split a movestack
buffalo.split <- split(buffalo)
class(buffalo.split)
class(buffalo.split[[1]])
names(buffalo.split)

## stack a list of move objects 
buffalo.stk <- moveStack(buffalo.split, forceTz="UTC")
buffalo.stk@timestamps[1]

# if argument forceTz is not stated, the timestamp is converted to the computer timezone


########################
#### OUTPUTTING DATA ###
########################

## save the move object for later
save(buffalo, file="buffalo_cleaned.Rdata")

## save as a text file
buffaloDF <- as.data.frame(buffalo)
write.table(buffaloDF, file="buffalo_cleaned.csv", sep=",", row.names = FALSE)

## save as a shape file
writeOGR(buffalo, getwd(), layer="buffalo", driver="ESRI Shapefile")


## write kml or kmz of a Movestack ##
library("plotKML")
# open a file to write the content
kml_open('buf.kml')
# write the movement data individual-wise
for(i in levels(trackId(buffalo)))
  kml_layer(as(buffalo[[i]],'SpatialLines'))
# close the file
kml_close('buf.kml')

## and export KML using writeOGR ##
for(id in levels(trackId(buffalo))){
  writeOGR(as(buffalo[[id]], "SpatialPointsDataFrame"),
           paste(id, ".kml", sep=""),
           id, driver="KML")
  
  writeOGR(as(buffalo[[id]], "SpatialLinesDataFrame"),
           paste(id, "-track.kml", sep=""),
           id, driver="KML")
  
  print(paste("Exported ", id,
              " successfully.", sep=""))
}


#########################
### NON LOCATION DATA ###
#########################

## Download non-location data as a data.frame 
acc <- getMovebankNonLocationData(study="MPIAB white stork lifetime tracking data (2013-2014)",
                                  sensorID="Acceleration",
                                  animalName="DER AR439", login=cred)
str(acc)

## Check sensors available in a specific study
getMovebankSensors(study="MPIAB white stork lifetime tracking data (2013-2014)", login=cred)[1:10,]

## List of all available sensor types on Movebank
getMovebankSensors(login=cred)[,3:5]

## visualize and get basic stats of acceleration data (currently only for eObs tags)
## moveACC : https://gitlab.com/anneks/moveACC


## NOTE: for more details see vignettes https://bartk.gitlab.io/move/articles/move.html
