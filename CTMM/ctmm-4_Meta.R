
###########################################################################
# Meta-analyses of population-level mean parameters #######################
###########################################################################

# Load the pre-requisite package:
library(ctmm)

# Load the buffalo tracking dataset:
data("buffalo")

# Load pre-run objects (if needed):
# load("meta.RData")

# fit movement models
FITS <- list()
for(i in 1:length(buffalo)) {
  GUESS <- ctmm.guess(buffalo[[i]], interactive = FALSE)
  FITS[[i]] <- ctmm.fit(buffalo[[i]], GUESS, trace = 2)
}

# calculate AKDEs on a consistent grid
AKDES <- list()
AKDES <- akde(buffalo, FITS, trace = 2)

# color to be spatially distinct
COL <- color(AKDES, by = "individual")

# plot AKDEs
plot(AKDES,
     col.DF = COL,
     col.level = COL,
     col.grid = NA,
     level = NA)

# cluster-analysis of buffalo
cluster(AKDES, sort = TRUE)

# meta-analysis of buffalo home-range areas
## What is the mean home range area of an average individual:

meta(AKDES,
     col = c(COL,"black"), 
     verbose = TRUE, # verbose output with CIs
     sort = TRUE) 

## Compare groups
meta(list(south = AKDES[1:3],
          north = AKDES[4:6]),
     plot = TRUE, 
     verbose = TRUE) 


# Population models / population range: -----------------------------------

MEAN.FITS <- mean(FITS)
summary(MEAN.FITS) 

## What is the population range?
MEAN <- mean(AKDES) # distribution of the sample
plot(buffalo, MEAN)

PKDE <- pkde(buffalo, AKDES) # distribution of the population
plot(buffalo, PKDE)

EXT <- extent(list(MEAN, PKDE))
COL <- c("red", "black", "blue")

par(mfrow = c(1,2))
plot(buffalo, MEAN, col = COL, ext = EXT)
title("mean()")
plot(buffalo, PKDE, col = COL, ext = EXT)
title("pkde()")
par(mfrow = c(1,1))

summary(MEAN)$CI
summary(PKDE)$CI 

# save(FITS = FITS,
#      AKDES = AKDES,
#      PKDE = PKDE,
#      file = "meta.RData")
