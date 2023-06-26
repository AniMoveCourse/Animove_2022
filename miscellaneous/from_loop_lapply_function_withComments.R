# tricks and tips:
# ALT+ - produces <- 
# CTRL+SHIFT+C produces #, you can also select an entire chunk of code an comment or uncomment it this way
# CTRL+i on a selected code, produces intention


# how to go from a script, to a for loop, to a lapply, to a function
library(move)
library(lubridate)
data(fishers)
# we make a calculation that works on one instance, e.g one individual of a movestack
myspeed <- speed(fishers[['Leroy']])
TL <- timeLag(fishers[['Leroy']], "mins")
dist <- distance(fishers[['Leroy']])
DF <- data.frame(date=date(fishers[['Leroy']]@timestamps[-1]), speed=myspeed, timeLag=TL, steplength=dist)
meanDF_leroy <- data.frame(aggregate(DF[ ,c('speed','timeLag','steplength')], by=list(DF$date), FUN=mean))
colnames(meanDF_leroy) <- c('Date', 'avgSpeed', 'avgTimelag', 'avgStepLength')

######################
## make a for loop 
########################
# now we want to use the upper piece of code to run through all our instances, e.g. individuals, so we can do a for loop to loop around each individual
mynames <- namesIndiv(fishers)
for(x in mynames){
  myspeed <- speed(fishers[[x]])
  TL <- timeLag(fishers[[x]], "mins")
  dist <- distance(fishers[[x]])
  DF <- data.frame(date=date(fishers[[x]]@timestamps[-1]), speed=myspeed,timeLag=TL,steplength=dist)
  meanDF <- data.frame(aggregate(DF[,c('speed','timeLag','steplength')], by=list(DF$date), FUN=mean))
  colnames(meanDF) <- c('Date','avgSpeed','avgTimelag','avgStepLength')
  meanDF$indiv <- x # we add this line, because we want to know which data belongs to each individual
  if(x==mynames[1]){meanDFall_loop <- meanDF}else{meanDFall_loop <- rbind(meanDFall_loop,meanDF)} # getting tables of each round rbinded, to in the end obtain a large table with all individuals
}

#####################
## make a lapply 
#####################
# a lapply is much more efficient than a for loop. In most cases a for loop can be easily transformed into a lapply. but there are cases where one just needs a for loop, and there is no way, or only a hard one, around it
mynames <- namesIndiv(fishers)
meanL <- lapply(mynames, function(x){
  myspeed <- speed(fishers[[x]])
  TL <- timeLag(fishers[[x]], "mins")
  dist <- distance(fishers[[x]])
  DF <- data.frame(date=date(fishers[[x]]@timestamps[-1]), speed=myspeed,timeLag=TL,steplength=dist)
  meanDF <- data.frame(aggregate(DF[,c('speed','timeLag','steplength')], by=list(DF$date), FUN=mean))
  colnames(meanDF) <- c('Date','avgSpeed','avgTimelag','avgStepLength')
  meanDF$indiv <- x
  return(meanDF) # you always need to tell a function what it should export as a result. you can only export ONE single object of what ever class and size
})
meanDFall_lapply <- do.call('rbind',meanL) # here we join the tables contained in the list into one big table

##########################################
### make a lapply using a list as an input
##########################################
# lapplys are designed to work with lists (there are other functions, apply, sapply, tapply, mapply that work on different data classes)
# so instead of a vector, we actually feed in a list, making the calculation more general, in this case it can be applied to any movestack
meanL <- lapply(split(fishers), function(x){
  myspeed <- speed(x)
  TL <- timeLag(x, "mins")
  dist <- distance(x)
  DF <- data.frame(date=date(x@timestamps[-1]), speed=myspeed,timeLag=TL,steplength=dist)
  meanDF <- data.frame(aggregate(DF[,c('speed','timeLag','steplength')], by=list(DF$date), FUN=mean))
  colnames(meanDF) <- c('Date','avgSpeed','avgTimelag','avgStepLength')
  meanDF$indiv <- namesIndiv(x)
  return(meanDF)
})
meanDFall_lapply_L <- do.call('rbind',meanL)

######################################
## make a function 
###############################
# when we do a lapply we automatically are creating a function. but of course we can also use a function in many other situations
# here we just copy/pasted the code from the above lapply and gave it a name
calcMeansPerDay <- function(x){
  myspeed <- speed(x)
  TL <- timeLag(x, "mins")
  dist <- distance(x)
  DF <- data.frame(date=date(x@timestamps[-1]), speed=myspeed,timeLag=TL,steplength=dist)
  meanDF <- data.frame(aggregate(DF[,c('speed','timeLag','steplength')], by=list(DF$date), FUN=mean))
  colnames(meanDF) <- c('Date','avgSpeed','avgTimelag','avgStepLength')
  meanDF$indiv <- namesIndiv(x)
  return(meanDF)
}

############################################
## lapply with a externally coded function
##############################################
# here we do exactly the same as before, but using the function that we created above
meanL <- lapply(split(fishers), calcMeansPerDay)
meanDFAll_function <- do.call('rbind',meanL)

###################################################
## create a function with more than one variable 
##############################################
# here in this case we create a function where we can input any move object, and calculate the time lag in our units of choice. in the previous calculations it was always in minutes

calcMeansPerDay_v2 <- function(x, units){ # x= is a move object; units= sec, mins, hours, days, etc
  myspeed <- speed(x)
  TL <- timeLag(x, units)
  dist <- distance(x)
  DF <- data.frame(date=date(x@timestamps[-1]), speed=myspeed,timeLag=TL,steplength=dist)
  meanDF <- data.frame(aggregate(DF[,c('speed','timeLag','steplength')], by=list(DF$date), FUN=mean))
  colnames(meanDF) <- c('Date','avgSpeed','avgTimelag','avgStepLength')
  meanDF$indiv <- namesIndiv(x)
  return(meanDF)
}

# here we are using the function created above
calcMeansPerDay_v2(x=fishers[["Leroy"]], units="hours")
calcMeansPerDay_v2(x=fishers[["Ricky.T"]], units="secs")
# or
meanL_v2 <- lapply(split(fishers), calcMeansPerDay_v2, units="hours")
meanDFAll_function_v2 <- do.call('rbind',meanL_v2)

## nesting several functions
do.call("rbind", lapply(split(fishers), calcMeansPerDay_v2, units="hours"))


## also have a look at the library(foreach) and library(doMC) and library(plyr)


## example to run code parallel across cores (this is just one way of doing it) 
library(doParallel)
library(plyr)
mycores <- detectCores()-1 # always make sure to not use all cores of your PC, at least leave 1 free
registerDoParallel(mycores) 

meanL_parallel <- llply(split(fishers), # your list
                        calcMeansPerDay, #your function (and parameters if necessary)
                        .parallel = T)


