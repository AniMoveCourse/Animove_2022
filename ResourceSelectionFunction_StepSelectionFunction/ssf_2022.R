#' # Step selection functions 
#'  Björn Reineking, 2022-09-13
#'  
#'  #' ## Preamble
#' - Movement is about solving the problem of life
#' - Trajectory analysis is about finding the problem
#'  
#' ## Step selection functions: Literature
#' - Forester et al. (2009) Ecology 90: 3554–3565
#' - Potts et al. (2014) Ecology and Evolution 4: 4578–4588
#' - Potts et al. (2014) J. R. Soc. Interface 11: 20140333
#' - Avgar et al. (2016) Methods in Ecology and Evolution 7: 619–630
#' 
#' ## Step selection functions: A receipe
#' 
#'- Problem formulation
#'- Data collection
#'- Data pre-processing
#'- Modelling
#'- Interpretation
#'- Iterate process
#'
#' We will use the buffalo data set.
#'
#' ## Loading packages
#+ results='hide', message=FALSE, warning=FALSE
library(ggplot2)
library(RStoolbox)
library(animove)
library(survival)
library(MASS)
library(lubridate)
library(dplyr)
library(nlme)
library(pbs)
library(circular)
library(CircStats)
library(amt)
library(ctmm)
library(move)

#' ## Load buffalo data
#' See Kami's slides for how we got here...
data(buffalo_utm)

#' ## Convert MoveStack to amt::track object
#' Currently, there is no conversion function from move::moveStack to amt::track implemented, so we do it by hand
buffalo_tracks <- as.data.frame(buffalo_utm) %>% 
  select(id = individual.local.identifier, x = coords.x1, y = coords.x2, ts = timestamp) %>% 
  make_track(x, y, ts, id, crs = projection(buffalo_utm))

#' ## Always inspect your data: summary statistics
summary(buffalo_tracks)


#' ## Environmental data: topography, waterways, and NDVI
data(buffalo_env)

#' We can plot the buffalo tracks, but they do not have a plot method, so 
#' we need to give a bit more information to the points() function.
raster::plot(raster(buffalo_env, 1))
points(y_ ~ x_, data = buffalo_tracks, col = factor(buffalo_tracks$id))

#' To speed up analyses, we will only work with one individual, Cilla
cilla <- filter(buffalo_tracks, id == "Cilla")

#' ## Inspect data
hist(step_lengths(cilla))
which(step_lengths(cilla)>5000)

#' The very first step is unusually long; let us plot the first day in red on top of the full trajectory.
plot(cilla, type = "l")
lines(filter(cilla, t_ < min(t_) + days(1)), col = "red")

#' Let us exclude the full first day
cilla <- filter(cilla, t_ > min(t_) + days(1))

#' ## Thin movement data and split to bursts
#' - We reduce the data set to observations that are within a certain time step range. The SSF assumes Brownian motion, so we should thin sufficiently, so that the velocities of successive steps are uncorrelated. See presentation by Chris on Monday. Here we go for 3 hours. 
#' - There is some tolerance around the target time interval of 3 hours. When two observations are separated by less than the threshold, the second observation is removed
#' - When two observations are separated by more than the upper threshold, the observations are assigned to different bursts.
#' 
#' It is a good idea to perform the analysis at several temporal scales, i.e. different step durations.
#' 
#' The initial sampling rate of cilla is about 1 hour:
#' 
summarize_sampling_rate(cilla)

#' ## Prior analysis of spatio-temporal autocorrelation
#' SSF assumes that the velocities are not autocorrelated. So we cannot simply do the analysis
#' at the highest temporal resolution of the data.
#' We could try to use the downweighting trick that we use for the RSF, but for the SSF, 
#' the time between successive steps will also affect the point estimates of the parameters,
#' so the problem is a bit more tricky.
#' As a rule-of-thumb, I would downsample the data such that the autocorrelation in velocities
#' has decayed to something between 1% and 2%.
#' Say you had a sampling of 1 hour, and the SSF should be done at 3 hours given this
#' rule of thumb, then you could still use all data, by constructing 3 step data sets, each starting one hour apart, and
#' give each a weight of 1/3 in the likelihood. 
#' 
#' 
library(ctmm)
cilla_telemetry <- as_telemetry(cilla)
plot(variogram(cilla_telemetry), xlim = c(0, 10 * 3600))

GUESS <- ctmm.guess(cilla_telemetry, interactive=FALSE)
FIT <- ctmm.fit(cilla_telemetry, GUESS)

plot(variogram(cilla_telemetry), xlim = c(0, 10 * 3600))
abline(v = FIT$tau["velocity"]  * -log(0.01) / 3600, col = "blue")
abline(v = FIT$tau["velocity"]  * -log(0.02) / 3600, col = "red")
legend("bottomright", lty = 1, col = c("blue", "red"), legend = c("1%", "2%"), title = "Velocity\nautocorrelation", 
       bty = "n")

#' Now we resample to 3 hour intervals, with a tolerance of 15 minutes
step_duration <- 3
cilla <- track_resample(cilla, hours(step_duration), tolerance = minutes(15))

#' Look at the new sampling rate
summarize_sampling_rate(cilla)

#' So there is at least two observations that are more than 3 hours 15 minutes apart, so there should be at least two bursts: 
table(cilla$burst_)

#' If there are bursts, we may want to filter bursts with very few locations. For example, to calculate a turning angle, we need at least three locations. So we often will want to filter out bursts with at least 3 observations:
cilla <- filter_min_n_burst(cilla, 3)

#' Convert locations to steps. We will have fewer rows in the step data frame than in the track data frame because the final position is not a complete step.
ssf_cilla <- steps_by_burst(cilla)

#' We still have steps without a turning angle (the first step in a burst)
which(is.na(ssf_cilla$ta_))
ssf_cilla <- filter(ssf_cilla, !is.na(ta_))

#' ## Empirical distances and turning angles
par(mfrow = c(1, 2))
hist(ssf_cilla$sl_, breaks = 20, main = "", 
  xlab = "Distance (m)")
hist(ssf_cilla$ta_,  main="",breaks = seq(-pi, pi, len=11),
      xlab="Relative angle (radians)")

#' ## Fit gamma distribution to distances
fexp <- fitdistr(ssf_cilla$sl_, "exponential")
fgamma <- amt::fit_distr(ssf_cilla$sl_, "gamma")
par(mfrow = c(1, 1))
hist(ssf_cilla$sl_, breaks = 50, prob = TRUE, 
     xlim = c(0, 8000), ylim = c(0, 2e-3),
     xlab = "Step length (m)", main = "")
plot(function(x) dexp(x, rate = fexp$estimate), add = TRUE, from = 0.1, to = 8000, col = "red")
plot(function(x) dgamma(x, shape = fgamma$params$shape,
                        scale = fgamma$params$scale), add = TRUE, from = 0.1, to = 8000, col = "blue")
legend("topright", col = c("red", "blue"), lty = 1,
       legend = c("exponential", "gamma"), bty = "n")

#' ## Fit von Mises distribution to angles
fvmises <- fit_distr(ssf_cilla$ta_, "vonmises")
par(mfrow = c(1, 1))
hist(ssf_cilla$ta_, breaks = 50, prob = TRUE, 
     xlim = c(-pi, pi),
     xlab = "Turning angles (rad)", main = "")
plot(function(x) dvonmises(x, mu = 0, kappa = fvmises$params$kappa), add = TRUE, from = -pi, to = pi, col = "red")


#' Create random steps. We typically get a warning that "Step-lengths or turning angles contained NA, which were removed", because of the missing turning angles at the start of a burst.
set.seed(2)
ssf_cilla <- steps_by_burst(cilla)
ssf_cilla <- random_steps(ssf_cilla, n_control = 200)

#' ## Sanity check: plot the choice set for a given step
my_step_id <- 3
ggplot(data = filter(ssf_cilla, step_id_ == my_step_id | (step_id_ %in% c(my_step_id - 1, my_step_id - 2) & case_ == 1)),
       aes(x = x2_, y = y2_)) + geom_point(aes(color = factor(step_id_))) + geom_point(data = filter(ssf_cilla, step_id_ %in% c(my_step_id, my_step_id - 1, my_step_id - 2) & case_ == 1), aes(x = x2_, y = y2_, color = factor(step_id_), size = 2))

#' ## Extract environmental covariates
#' I recommend to always use the option "both", which provides the environmental conditions at the start and the end of the step.
#' The condition at the end are what we use for selection, and the conditions at the start can be used
#' to modify e.g. turning angles and step length.
ssf_cilla <- extract_covariates(ssf_cilla, buffalo_env, where = "both")

#' ## Add variable hour
#' Adding hour modelling diurnal variation in step lengths, turning angles, and preference for environmental conditions
ssf_cilla <- mutate(ssf_cilla, "hour" = hour(t1_) + minute(t1_) / 60)

#' Remove NA's
ssf_cilla <- ssf_cilla[complete.cases(ssf_cilla),]

#' ## A first model
m_1 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + slope_end + elev_end + 
                    mean_NDVI_end + var_NDVI_end + strata(step_id_))
summary(m_1)

#' ## Collinearity
#' In statistics, multicollinearity (also collinearity) is a phenomenon in which one predictor variables in a multiple regression model can be linearly predicted from the others with a substantial degree of accuracy. In this situation the coefficient estimates of the multiple regression may change erratically in response to small changes in the model or the data." [Wikipedia, accessed 29.08.2017](https://en.wikipedia.org/wiki/Multicollinearity)
#' 
#' One way of dealing with collinearity is to select a subset of variables that is sufficiently uncorrelated [Dormann et al. 2013](http://onlinelibrary.wiley.com/doi/10.1111/j.1600-0587.2012.07348.x/abstract). Here we simply look at pairwise correlation between predictors.
#' 
round(cor(ssf_cilla[, c("slope_end", "elev_end", 
  "water_dist_end", "mean_NDVI_end", "var_NDVI_end")]), 2)

#' elev and water_dist are positively correlated > 0.7

#' ## Which collinear variable to pick? 
#' - The one that is more relevant
#' - The one that is by itself a better predictor
#' 

m1_water <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + water_dist_end + strata(step_id_))
m1_elev <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + strata(step_id_))
AIC(m1_water$model)
AIC(m1_elev$model)

#' So we pick elev, because it by itself explains the movement better

#' Fit step selection function
m_1 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + slope_end + elev_end + 
                    mean_NDVI_end + var_NDVI_end + strata(step_id_))
summary(m_1)

#' slope and var_NDVI do not contribute significantly to the fit

#' ## Model selection
#' Model selection is a vast topic. I recommend using only few models with ecological justification, rather than 
#' searching for the "best" model in a huge model space.
#' Here we just use stepwise backward selection based on AIC

m_2 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + slope_end + elev_end + mean_NDVI_end + strata(step_id_))
AIC(m_1$model)
AIC(m_2$model)
summary(m_2)

m_3 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + mean_NDVI_end + strata(step_id_))
AIC(m_3$model)
summary(m_3)

#' ## Model checking: serial autocorrelation
#' Forester et al. 2009 Ecology 90:3554–3565.
#' Calculate  deviance residuals for each stratum (i.e., the sum of the residuals for the case and all associated controls).
#' 
ssf_residuals <- function(m, data) {
  df <- tibble::as_tibble(data.frame("time" = data$t1_,
                                       "residuals" = residuals(m$model, type = "deviance")))
  df <- df %>% dplyr::group_by(time) %>% dplyr::summarise(residuals = sum(residuals))
  df$group <- 1
  df
}

resid_df <- ssf_residuals(m_3, ssf_cilla)

#' Fit an intercept-only mixed-effects model using lme() from the nlme package.
#' 
rm1 <- lme(residuals ~ 1, random = ~ 1 | group, 
           data = resid_df)
plot(ACF(rm1, maxLag = 40), alpha = 0.05)

#' So we see that there is some significant autocorrelation at lag 1.
#' 
#' One effect of residual temporal autocorrelation is too extreme p-values, but it may also cause bias in parameter estimates.
#' Forester et al. 2009 suggest a way to estimate more appropriate confidence intervals and p-values.
 
#' ## Model evaluation
#' - R2 is low. Always is.
#' - Not yet clear what a good performance index would be. Perhaps https://github.com/aaarchmiller/uhcplots will help.
#' - Cross-validation
#'     - split steps in e.g. 5 folds (long stretches better - should be long enough so that autocorrelation in residuals has tapered off)
#'     - leave each fold out, refit and predict to left-out fold
#'     


#' ## Interpretation
#'  - Map preference
#' - Response functions
#'

#'## Map habitat preference ("Habitat map")
#' Caveat: Note that this habitat preference in general will not match the utilisation distribution of the animal (i.e. how much time it spends where). See below for more information.
#' The raster prediction function assumes that all environmental layers are
#' represented in one raster stack
#'
#' We need to exclude the parameters that are not related to the environment, i.e. the first three parameters related to turning angles and step length
coef(m_1)
remove_end <- function(x) {
  names(x) <- gsub("_end", "", names(x))
  x
}
habitat_map <- habitat_kernel(remove_end(coef(m_1)[-(1:3)]), buffalo_env, exp = FALSE)

ggp <- ggR(habitat_map, geom_raster = TRUE) +
  scale_fill_gradient(low = "lightgray", high = "black")
ggp + geom_path(data = ssf_cilla, aes(x = x1_, y = y1_))

#' Now zoom in
ggp <- ggR(crop(habitat_map, extent(as_sp(cilla)) + 5000), geom_raster = TRUE) +
  scale_fill_gradient(low = "lightgray", high = "black")
ggp + geom_path(data = ssf_cilla, aes(x = x1_, y = y1_, color = factor(burst_)))

#' The model is strongly driven by elevation
ggp <- ggR(crop(buffalo_env, extent(as_sp(cilla)) + 5000), layer = "elev", geom_raster = TRUE) +
  scale_fill_gradient(low = "lightgray", high = "black")
ggp + geom_path(data = ssf_cilla, aes(x = x1_, y = y1_, color = factor(burst_)))

#' ## Iterate
#' Here: a model with time-varying preference for mean_NDVI
#' We group observation in 3 hour bins to smooth the picture
#' 
boxplot(mean_NDVI_end ~ I(floor(hour/3)*3), data = 
  filter(ssf_cilla, case_ == 1),xlab = "Time of day", 
  ylab = "mean NDVI")

#' We can do the same for other variables, including those of the "movement kernel", e.g. distance
boxplot(sl_ ~ floor(hour), data = 
          filter(ssf_cilla, case_ == 1), xlab = "Time of day", 
        ylab = "dist")

#' What behavioural rhythm do these figures suggest?
m_time_ndvi <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + mean_NDVI_end + 
                   mean_NDVI_end:pbs(hour, df = 5, Boundary.knots = c(0,24)) + strata(step_id_))
m_time_dist <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + mean_NDVI_end + 
                   sl_:pbs(hour, df = 5, Boundary.knots = c(0,24)) + strata(step_id_))

#' ## Predictions with the model: response function
#' 
extract_stratum <- function(object) {
  attr(object$terms, 'special')$strata[1]
}
stratum <- extract_stratum(m_time_ndvi$model)
pred_data_ndvi <- data.frame("step_id_" = stratum, ta_ = 0, sl_ = 1, elev_end = 0, mean_NDVI_end = 1, hour = seq(0, 24, len = 101))
m_time_ndvi <- clogit(case_ ~ cos(ta_) + sl_ + log(sl_) + elev_end + mean_NDVI_end + 
                            mean_NDVI_end:pbs(hour, df = 5, Boundary.knots = c(0,24)) + strata(step_id_), data = ssf_cilla)
pred_time <- survival:::predict.coxph(m_time_ndvi, newdata = pred_data_ndvi, se.fit = TRUE)
upper <- pred_time$fit + 1.96 * pred_time$se.fit
lower <- pred_time$fit - 1.96 * pred_time$se.fit

par(mfrow = c(1, 1))
plot(pred_data_ndvi$hour, pred_time$fit, type = "l", 
  ylim = range(c(upper, lower)), xlab = "Time of day",
  ylab = "Preference mean_NDVI")
lines(pred_data_ndvi$hour, upper, lty = 2)
lines(pred_data_ndvi$hour, lower, lty = 2)
abline(h = 0, lty = 3)

#' ## Simulating with the model
set.seed(2)
k1 <- amt:::redistribution_kernel(m_3, map = buffalo_env, start = amt:::make_start.steps_xyt(ssf_cilla[1, ]),
                                  stochastic = TRUE, tolerance.outside = 0.2, as.raster = FALSE, 
                                  n.control = 1e3)
s1 <- amt:::simulate_path.redistribution_kernel(k1, n.steps = 500)

extent_tracks <- function(x, y) {
  df <- data.frame(na.omit(rbind(x[,c("x_", "y_")], y[,c("x_", "y_")])))
  raster::extent(c(range(df$x_), range(df$y_)))
}

elev_crop <- crop(buffalo_env[["elev"]], extent_tracks(s1, cilla) + 2000)
raster::plot(elev_crop)
lines(cilla)
lines(s1$x_, s1$y_, col = "red")

#' # Home ranging behaviour (Thanks, Chris!)
m_hr <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + 
                    elev_end + water_dist_end + x2_ + y2_ + I(x2_^2 + y2_^2) + strata(step_id_))
summary(m_hr)

k2 <- amt:::redistribution_kernel(m_hr, map = buffalo_env, start = amt:::make_start.steps_xyt(ssf_cilla[1, ]),
                            stochastic = TRUE, tolerance.outside = 0.01, as.raster = FALSE,
                            n.control = 1e3)
set.seed(2)
s2 <- amt:::simulate_path.redistribution_kernel(k2, n.steps = 1000)

raster::plot(crop(buffalo_env[["elev"]], extent_tracks(s2, cilla) + 2000))
lines(cilla)
lines(s2$x_, s2$y_, col = "red")

#' ## Barriers
water <- buffalo_env[["water_dist"]] > 100
water <- crop(water, extent(water) - 5000)
raster::plot(water)
ww <- clump(water)
raster::plot(ww)
ww[ww == 0] <- 2
names(ww) <- "water_crossed"
buffalo_env_2 <- raster::stack(ww, crop(buffalo_env, ww))

set.seed(2)
ssf_cilla <- cilla %>% steps_by_burst() %>% random_steps(n_control = 1000) %>% 
  extract_covariates(buffalo_env_2, where = "both") %>% 
  mutate(hour = hour(t1_) + minute(t1_) / 60) %>% 
  filter(complete.cases(.))

m_crossing <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + 
                    elev_end + water_dist_end + x2_ + y2_ + I(x2_^2 + y2_^2) +
                    I(water_crossed_end != water_crossed_start) + strata(step_id_))

summary(m_crossing)

set.seed(2)
k3 <-  amt:::redistribution_kernel(m_crossing, map = stack(buffalo_env_2), 
                                   start = amt:::make_start.steps_xyt(ssf_cilla[1, ]),
                            stochastic = TRUE, tolerance.outside = 0.01, 
                            as.raster = FALSE, n.control = 1e3)

s3 <- amt:::simulate_path.redistribution_kernel(k3, n.steps = 1000)

raster::plot(crop(buffalo_env[["elev"]], extent_tracks(s3, cilla) + 2000))
lines(cilla)
lines(s3$x_, s3$y_, col = "red")


#' ## From preference to utilisation maps
#' There is a fast method if we have a symmetric and temporally stable jump kernel, e.g. exponential, and no effect of step angles:
#' Barnett, A. & Moorcroft, P. (2008) Analytic steady-state space use patterns and rapid computations in mechanistic home range analysis. [Journal of Mathematical Biology, 57, 139–159](https://link.springer.com/article/10.1007/s00285-007-0149-8).
#' The generic but computationally expensive method is to do simulations: Signer et al. (2017) Estimating utilization distributions from fitted step-selection functions. [Ecosphere 8: e01771](http://onlinelibrary.wiley.com/store/10.1002/ecs2.1771/asset/ecs21771.pdf?v=1&t=j6xal7ze&s=ad768244a70741207c1e663409e7ab7e201c7132)

#' # Dependence of results on step interval
#' 
step_durations <- 3:12
do_run <- TRUE
buffalo_ids <- levels(factor(buffalo_tracks$id))# c("Cilla", "Gabs", "Mvubu")# levels(factor(buffalo_tracks$id))
if(do_run) {
step_interval_simulation_amt <- lapply(buffalo_ids, function(animal) {
  lapply(step_durations, function(step_duration) {
    ssf_animal <- filter(buffalo_tracks, id == animal) %>% 
      track_resample(hours(step_duration), tolerance = minutes(15)) %>%
      filter_min_n_burst(3) %>%
      steps_by_burst() %>% 
      random_steps(n = 200) %>%
      extract_covariates(buffalo_env) %>%
      filter(complete.cases(.)) %>% filter(sl_ > 0)
    m_1 <- clogit(case_ ~ cos(ta_) + sl_ + log(sl_) + elev + 
                    mean_NDVI + strata(step_id_), data = ssf_animal)
    list("coef" = coef(m_1), "confint" = confint(m_1))
  })
})
}

#' ## Parameter estimates
do_run <- TRUE
if (do_run)  {
for (j in seq(step_interval_simulation_amt)) {
  model_list <- step_interval_simulation_amt[[j]]
  name <- names(split(buffalo_utm))[j]
  coefs <- sapply(model_list, function(x) x$coef)
  ci_lower <- sapply(model_list, function(x) x$confint[, 1])
  ci_upper <- sapply(model_list, function(x) x$confint[, 2])
  par(mfrow = c(3, 2), mar=c(4,5,1,1), oma = c(1,1,5,1))
  for (i in rownames(coefs)) {
    plot(c(0,0), xlim = range(step_durations),
         ylim = range(c(0, ci_lower[i, ],
                        coefs[i,],
                        ci_upper[i, ])), type = "n", xlab = "Step duration (hr)", 
         ylab = i)
    abline(h = 0, lty = 2, col = "red")
    lines(step_durations, ci_lower[i, ], lty = 3)
    lines(step_durations, ci_upper[i, ], lty = 3)
    lines(step_durations, coefs[i, ])
  }
  mtext(name, outer = TRUE)
}
}

#' ## Some of the stuff to expand on
#' - Incorporating boundaries (e.g. rivers)
#' - Correct confidence intervals (Forrester et al. 2009)
#' - Efficient estimate of utilisation distribution
#' - Model performance
#' - Interactions between individuals
#' - Mixed effects models -> see new manuscript by Johannes Signer https://www.biorxiv.org/content/early/2018/09/08/411801
#' - Memory (or "future expectation")
#' - Interaction with changes of internal states (coming from another model/other observations)
#' 
