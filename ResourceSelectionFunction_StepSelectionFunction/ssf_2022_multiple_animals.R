#' # Step selection functions: multiple animals
#' Björn Reineking, 2022-09-14
#' 
#' We will use the buffalo data set.
#' The code in this example is based on (slightly modified) from Muff, Signer & Fieberg (2020) J Anim Ecol 89: 80-92.
#' 
#' # Loading packages
#+ results='hide', message=FALSE, warning=FALSE
library(animove)
library(dplyr)
library(amt)
library(move)
library(glmmTMB)

#' Load data
data(buffalo_utm)
#' Convert MoveStack to amt::track object.
#' 
#' Currently, there is no conversion function from move::moveStack to amt::track implemented, so we do it by hand
buffalo_tracks <- amt::mk_track(as.data.frame(buffalo_utm),
                        coords.x1, coords.x2, timestamps, 
                        id = individual.local.identifier,
                        crs = crs(buffalo_utm))
#' Environmental data: topography, waterways, and NDVI
data(buffalo_env)
#' We exclude the full first day, since at least for buffalo Cilla there are some weird locations
buffalo_filtered <- buffalo_tracks %>% group_by(id) %>% filter(t_ > min(t_) + days(1))
#' ## Thin movement data and split to bursts
#' 
#' - We reduce the data set to observations that are within a certain time step range. The SSF assumes Brownian motion, so we should thin sufficiently, so that the velocities of successive steps are uncorrelated. Here we go for 3 hours. 
#' - There is some tolerance around the target time interval of 3 hours. When two observations are separated by less than the threshold, the second observation is removed
#' - When two observations are separated by more than the upper threshold, the observations are assigned to different bursts.
#' 
step_duration <- 3
buffalos <- track_resample(buffalo_filtered, hours(step_duration), tolerance = minutes(15))

#' Often will want to filter out bursts with at least 3 observations:
buffalos <- filter_min_n_burst(buffalos, 3)
#' Convert locations to steps. We will have fewer rows in the step data frame than in the track data frame because the final position is not a complete step. 
buffalos <- ungroup(buffalos)
buffalo_list <- split(buffalos, buffalos$id)
ssf_buffalos <- map(buffalo_list, steps_by_burst)
ssf_buffalos <- bind_rows(ssf_buffalos, .id = "id")

#' Create random steps. We typically get a warning that "Step-lengths or turning angles contained NA, which were removed", because of the missing turning angles at the start of a burst.

set.seed(2)
ssf_buffalos <- random_steps(ssf_buffalos, n_control = 100)

#' Extract environmental covariates
ssf_buffalos <- extract_covariates(ssf_buffalos, buffalo_env)
#' Remove NA's
ssf_buffalos <- ssf_buffalos[complete.cases(ssf_buffalos),]
#' Convert to standard data frame
ssf_df <- as.data.frame(ssf_buffalos)

#' ## Model without taking individal variation into account
m_3 <- fit_clogit(ssf_df, case_ ~ cos(ta_) + sl_ + log(sl_) + 
                    elev + mean_NDVI + strata(step_id_))
summary(m_3)

#' ## Fit separate parameter values for each individual
m_4 <- fit_clogit(ssf_df, case_ ~ (cos(ta_) + sl_ + log(sl_) + 
                                     elev + mean_NDVI):id + strata(step_id_))
summary(m_4)

#' ## Fixed-effects model glmmTMB
#' To see if we get (more or less) the same result with the glmmTMB library
#' 
TMBStruc.fix = glmmTMB(case_ ~ -1 + cos(ta_) + sl_ + log(sl_) + 
                         elev + mean_NDVI  + (1|step_id_), 
                       family=poisson, data=ssf_df, doFit=FALSE) 
#' Fix the standard deviation of the first random term, which is the `(1|step_id_)` component  in the above model equation:
TMBStruc.fix$parameters$theta[1] = log(1e3) 

#' We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
TMBStruc.fix$mapArg = list(theta=factor(c(NA)))

#' Then fit the model and look at the results:
glmm.TMB.fixed <- glmmTMB:::fitTMB(TMBStruc.fix)

#' OK, the model did not converge. So we have to tweak the options for the optimizer:
#' 

TMBStruc.fix_2 = glmmTMB(case_ ~ -1 + cos(ta_) + sl_ + log(sl_) + 
                         elev + mean_NDVI  + (1|step_id_), 
                       family=poisson, data=ssf_df, doFit=FALSE, 
                       control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))) 
TMBStruc.fix_2$parameters$theta[1] = log(1e3) 
TMBStruc.fix_2$mapArg = list(theta=factor(c(NA)))
glmm.TMB.fixed_2 <- glmmTMB:::fitTMB(TMBStruc.fix_2)
summary(glmm.TMB.fixed_2)
confint(glmm.TMB.fixed_2)

#' ## Random effects model
#' 
#' The variable representing different animals, ANIMAL_ID, needed to be numeric in earlier versions of glmmTMB;
ssf_df$ANIMAL_ID <- factor(ssf_df$id)

TMBRandomStruc = glmmTMB(case_ ~ -1 + cos(ta_) + sl_ + log(sl_) + elev + mean_NDVI + 
                           (1|step_id_) + 
                           (0 + cos(ta_) | ANIMAL_ID) + 
                           (0 + sl_ | ANIMAL_ID) + 
                           (0 + log(sl_) | ANIMAL_ID) + 
                           (0 + elev | ANIMAL_ID) +
                           (0 + mean_NDVI | ANIMAL_ID), family = poisson, 
                         data=ssf_df, doFit=FALSE, 
                         control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

TMBRandomStruc$parameters$theta[1] = log(1e3) 

#' Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
TMBRandomStruc$mapArg = list(theta=factor(c(NA,1:5)))

#' Fit the model and look at the summary:
glmm.TMB.random <- glmmTMB:::fitTMB(TMBRandomStruc)

summary(glmm.TMB.random)
#' So elevation is significant at the population level (animals avoid higher elevation sites).
#' 

#' 95\% CIs for fixed and random effects (standard deviations) are obtained via the confint() function:
confint(glmm.TMB.random)
#' The standard deviation of the random effect is large relative to the fixed effect size,
#' so there are animals that react positively to elevation (as we had also seen for Toni in the model m_4 with separate parameter values for each individual).