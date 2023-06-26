
###########################################################################
# Autocorrelated home range estimation ####################################
###########################################################################

# Load the pre-requisite package:
library(ctmm)

# Load pre-run objects (if needed):
# load("akde.RData")

# Data --------------------------------------------------------------------

# Loading tracking datasets:
data("buffalo")
data("gazelle")

animal1_buffalo <- buffalo$Pepper # or buffalo[[4]]
head(animal1_buffalo)

animal2_gazelle <- gazelle[[11]]
head(animal2_gazelle)

# Plotting locations:
plot(animal1_buffalo, col = "red", lwd = 3)
plot(animal2_gazelle, col = "blue", lwd = 3)


# Range residency assumption ----------------------------------------------

level <- 0.95 # we want to display 95% confidence intervals
xlim <- c(0,1 %#% "day") # to create a window of one day

# Checking for the range residency assumption:
SVF <- variogram(animal1_buffalo)
par(mfrow = c(1,2))
plot(SVF, fraction = 0.5, level = level)
abline(v = 1, col = "red", lty = 2) # adding a line at 1 month
plot(SVF, xlim = xlim, level = level)
par(mfrow = c(1,1))


# Model selection ---------------------------------------------------------
# Selecting the best-fit movement model through model selection:

# Calculate an automated model guesstimate:
GUESS1 <- ctmm.guess(animal1_buffalo, interactive = FALSE)

# Automated model selection, starting from GUESS1:
start_time <- Sys.time()
FIT1_ML <- ctmm.select(animal1_buffalo, GUESS1, 
                       method = "ML", verbose = TRUE)
Sys.time() - start_time # Time difference of 2.435679 mins
summary(FIT1_ML)

plot(SVF, CTMM = FIT1_ML[[1]],
     units = TRUE, fraction = 0.5, level = c(0.95, 0.50), 
     col = "black", col.CTMM = "red")

start_time <- Sys.time()
FIT1_pHREML <- ctmm.select(animal1_buffalo, GUESS1,
                           method = "pHREML", verbose = TRUE)
## reminder: it will default to pHREML if no method is specified.
Sys.time() - start_time # Time difference of 2.435679 mins
summary(FIT1_pHREML)

plot(SVF, CTMM = FIT1_pHREML[[1]],
     units = TRUE, fraction = 0.5, level = c(0.95, 0.50), 
     col = "black", col.CTMM = "red")

summary(FIT1_ML[[1]]) # best-fit model only
summary(FIT1_pHREML[[1]])

plot(SVF, CTMM = FIT1_pHREML,
     units = TRUE, level = 0.95, 
     # xlim = c(0, 10 %#% "hours"),
     col = "black", 
     col.CTMM = rainbow(length(FIT1_pHREML)))

# Home range estimator ----------------------------------------------------
# Feeding a movement model into the home range estimator

# Run an area-corrected AKDE (default):
AKDE1_ML <- akde(animal1_buffalo, FIT1_ML, debias = TRUE)
AKDE1_pHREML <- akde(animal1_buffalo, FIT1_pHREML, debias = TRUE)

summary(AKDE1_ML, level.UD = 0.95)$CI # 95% home range area
summary(AKDE1_pHREML, level.UD = 0.95)$CI

( 1 - summary(AKDE1_ML)$CI[1,2] / summary(AKDE1_pHREML)$CI[1,2] ) * 100
# ML overestimates by 5%

m.iid <- ctmm.fit(animal1_buffalo) # IID
KDE1 <- akde(animal1_buffalo, m.iid)

# Creating an extent that includes both UDs at the 95% CI level:
newEXT <- extent(list(AKDE1_pHREML, KDE1))

# Plotting KDE and AKDE side-by-side:
par(mfrow = c(1,2))
plot(animal1_buffalo, UD = KDE1, ext = newEXT)
title(expression("KDEc"))
plot(animal1_buffalo, UD = AKDE1_pHREML, ext = newEXT)
title(expression("AKDEc"))
par(mfrow = c(1,1))


# Mitigation measures -----------------------------------------------------
# Evaluating additional biases, applying mitigation measures

## Irregular representation in time: --------------------------------------

plot(animal1_buffalo, lwd = 3)

## Sample sizes:
summary(AKDE1_pHREML)$DOF["area"] # effective sample size of animal1
nrow(animal1_buffalo) # absolute sample size

# plot all sampling intervals
dt.plot(animal1_buffalo) # Pepper (buffalo[[4]])
abline(h = 2 %#% "hours", col = "red")

dt.plot(buffalo$Cilla) # Cilla (buffalo[[1]])
summary(animal1_buffalo)

col <- "hr" %#% diff(animal1_buffalo$t)
# minimum adjacent sampling interval
col <- pmin(c(Inf,col),c(col,Inf))
# sampling intervals under 1.5 hours
col <- (col < 1.5)
# red (low-frequency) or yellow (high-frequency)
col <- grDevices::rgb(1, col, 0)
plot(animal1_buffalo, col = col, lwd = 2)

# minimum sampling interval
"minutes" %#% min(diff(animal1_buffalo$t))

# Option 1) exact calculation - slow on larger datasets
start_time <- Sys.time()
option1_wAKDE1 <- akde(animal1_buffalo,
                     CTMM = FIT1_pHREML,
                     weights = TRUE, 
                     fast = FALSE, PC = "direct")
Sys.time() - start_time # Time difference of 3.063514 mins

# Option 2) use tiny dt
dt <- min(diff(animal1_buffalo$t))
start_time <- Sys.time()
option2_wAKDE1 <- akde(animal1_buffalo,
                     CTMM = FIT1_pHREML,
                     weights = TRUE, dt = dt)
Sys.time() - start_time # Time difference of 43.70387 secs

# Option 3) remove tiny sampling intervals and use smaller dt (default)
dt <- 1 %#% "hr"
start_time <- Sys.time()
option3_wAKDE1 <- akde(animal1_buffalo,
                      CTMM = FIT1_pHREML,
                      weights = TRUE, dt = dt)
Sys.time() - start_time # Time difference of 1.384689 secs

summary(option1_wAKDE1)$CI
summary(option2_wAKDE1)$CI
summary(option3_wAKDE1)$CI

EXT <- extent(list(option1_wAKDE1, 
                   option2_wAKDE1, 
                   option3_wAKDE1), level = 0.95)
par(mfrow = c(1,3))
plot(animal1_buffalo, UD = option1_wAKDE1,
     col = col, ext = newEXT)
plot(animal1_buffalo, UD = option2_wAKDE1, 
     col = col, ext = newEXT)
plot(animal1_buffalo, UD = option3_wAKDE1, 
     col = col, ext = newEXT)
par(mfrow = c(1,1))

wAKDE1_pHREML <- option1_wAKDE1

summary(wAKDE1_pHREML)$CI # 95% home range area (weighted)

EXT <- extent(list(AKDE1_ML, AKDE1_pHREML, wAKDE1_pHREML), level = 0.95)

# Plotting pHREML (with and without weights) side-by-side:
par(mfrow = c(1,2))
plot(animal1_buffalo, UD = AKDE1_pHREML, ext = newEXT)
title(expression("pHREML AKDE"["C"]))
plot(animal1_buffalo, UD = wAKDE1_pHREML, ext = newEXT)
title(expression("pHREML wAKDE"["C"]))
par(mfrow = c(1,1))

( 1 - summary(AKDE1_pHREML)$CI[1,2] / summary(wAKDE1_pHREML)$CI[1,2] ) * 100


## Low sample sizes: ------------------------------------------------------

plot(animal2_gazelle, lwd = 3)

GUESS2 <- ctmm.guess(animal2_gazelle, interactive = FALSE)

FIT2_ML <- ctmm.select(animal2_gazelle, GUESS2, method = "ML")
FIT2_pHREML <- ctmm.select(animal2_gazelle, GUESS2, method = "pHREML")
summary(FIT2_pHREML)

# Expected order of ML bias:
1/summary(FIT2_ML)$DOF["area"]

UD2_ML <- akde(animal2_gazelle, FIT2_ML)
UD2_pHREML <- akde(animal2_gazelle, FIT2_pHREML)

summary(UD2_ML)$CI
summary(UD2_pHREML)$CI

( 1 - summary(UD2_ML)$CI[1,2] / summary(UD2_pHREML)$CI[1,2] ) * 100

summary(UD2_pHREML)$DOF["area"] # effective sample size
nrow(animal2_gazelle) # absolute sample size

# Expected order of pHREML bias:
1/summary(FIT2_pHREML)$DOF["area"]^2

start_time <- Sys.time() # start recording running time
BOOT <- ctmm.boot(animal2_gazelle, FIT2_pHREML, 
                  error = 0.01, trace = 2, cores = -1)
# save(BOOT, file = here::here("code", "outputs", "bootstrap.RData"))
## note: this function incurs substantial computational cost, may take hours.
( total_time <- Sys.time() - start_time ) # output running time
#Time difference of 43.93944 mins

summary(BOOT)
1/summary(BOOT)$DOF["area"]^3 # expected order of bias

UD2_bpHREML <- akde(animal2_gazelle, BOOT, weights = TRUE)

summary(UD2_pHREML)$CI
summary(UD2_bpHREML)$CI

( 1 - summary(UD2_pHREML)$CI[1,2] / summary(UD2_bpHREML)$CI[1,2] ) * 100

EXT <- extent(list(UD2_pHREML, UD2_bpHREML), level = 0.95)

# Plotting pHREML and bootstrapped-pHREML side-by-side:
par(mfrow = c(1,2))
plot(animal2_gazelle, UD = UD2_pHREML, ext = EXT)
title(expression("pHREML AKDE"["C"]))
plot(animal2_gazelle, UD = UD2_bpHREML, ext = EXT)
title(expression("Bootstrapped pHREML wAKDE"["C"]))

# save(FIT1_ML,
#      FIT1_pHREML,
#      wAKDE1_pHREML,
#      FIT2_ML,
#      FIT2_pHREML,
#      BOOT,
#      file = "akde.RData")
