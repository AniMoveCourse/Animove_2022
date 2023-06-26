
###########################################################################
# Home range overlap and encounter rates ##################################
###########################################################################

# Load the pre-requisite package:
library(ctmm)

# Load the buffalo tracking dataset:
data("buffalo")

# Load pre-run objects (if needed):
# load("overlap.RData")

# Subset the individuals of interest
ANIMALS <- buffalo[4:5]
plot(ANIMALS)

# And return some summary statistics
summary(ANIMALS)
head(ANIMALS$Pepper)
head(ANIMALS$Queen)

plot(ANIMALS, col = c("#000000", "#b30c00"))

# Estimate the empirical variograms
SVF_1 <- variogram(ANIMALS[[1]])
SVF_2 <- variogram(ANIMALS[[2]])

# Visually inspect these for range residency
par(mfrow = c(1,2))
plot(SVF_1)
plot(SVF_2)
par(mfrow = c(1,1))

# Generate initial guesses of the parameter estimates
GUESS_1 <- ctmm.guess(ANIMALS[[1]], interactive = FALSE)
GUESS_2 <- ctmm.guess(ANIMALS[[2]], interactive = FALSE)

#Fit and select movement models for each animal
FITS_1 <- ctmm.select(ANIMALS[[1]], GUESS_1)
FITS_2 <- ctmm.select(ANIMALS[[2]], GUESS_2)

# Return a summary of the selected model
summary(FITS_1) # for the first individual
summary(FITS_2) # for the second individual

# Store the best fit movement models in a list
FITS <- list(FITS_1, FITS_2)

# Estimate the HR estimates for each individual
HR_UDS <- akde(ANIMALS,
               FITS,
               weights = TRUE)

#Plot the data and HR estimates
plot(ANIMALS,
     UD = HR_UDS,
     col = c("#000000", "#b30c00"),
     col.DF = c("grey50", "#d60e00"),
     col.grid = NA)

# Estimate the home range overlap
OVERLAP <- overlap(HR_UDS, level = .95)
OVERLAP

# Estimate the encounter rate
CDE <- encounter(HR_UDS, level = .95)
summary(CDE)

plot(ANIMALS,
     UD = CDE,
     col = c("#000000", "#b30c00"),
     col.DF = "#009da0",
     col.grid = NA)
     
# save(ANIMALS = ANIMALS,
#      FITS_1 = FITS_1,
#      FITS_2 = FITS_2,
#      FITS = FITS,
#      HR_UDS = HR_UDS,
#      file = "overlap.RData")
