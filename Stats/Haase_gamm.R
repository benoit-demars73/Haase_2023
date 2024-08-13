# -------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphics windows 

# -------------------------------------------------------------------
#

# Noam Ross on spatio-temporal gam https://noamross.github.io/gams-in-r-course/chapter3
# Gavin Simpson on random effect https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/
# bs='re' to have a random effect
# Gavin Simpson on autocorrelation and gamm
# https://stats.stackexchange.com/questions/80823/do-autocorrelated-residual-patterns-remain-even-in-models-with-appropriate-corre

# GAMM for autocorrelation
# https://r.qcbs.ca/workshop08/book-en/introduction-to-generalized-additive-mixed-models-gamms.html

library(dplyr) 
library("lubridate")
library(mgcv)
library(gamm4)
library(ggplot2)

setwd("~/Haase_2023/Data")

# retrieve table produce in Taxonomic_richness_anomaly.R script
df.taxa.90 <- read.table("~/Haase_summary_taxa_rank_5y_1990-2020.txt", sep ="\t", header=TRUE)
head(df.taxa.90)
nrow(df.taxa.90) # 18657
n_distinct(df.taxa.90$site_id) # 1120
str(df.taxa.90)

df.taxa.90$year <- as.numeric(df.taxa.90$year) 
df.taxa.90$site_id <- as.factor(df.taxa.90$site_id) 
df.taxa.90$study_id <- as.factor(df.taxa.90$study_id) 
df.taxa.90$c.Fy_id_round <- as.factor(df.taxa.90$c.Fy_id_round) 

# another idea is just to provide corrected Species_richness anomaly (5y) as follows (see Fig.4 middle panel)
df.taxa.90$Species_richness_corr = df.taxa.90$Species_richness-(df.taxa.90$c.Fy_id*12.85)
df.taxa.90 <- transform(df.taxa.90, 
                  Species_richness_corr_5y = 
                  Species_richness_corr - 
                  ave(replace(Species_richness_corr, year < 2003 | year > 2007, NA), 
              site_id, FUN=function(t) mean(t, na.rm=TRUE)))
head(df.taxa.90)


# create iYear.90
df.taxa.90$iYear.90 <-  df.taxa.90$year - 1990 +1
df.taxa.90$iYear.90 <- as.factor(df.taxa.90$iYear.90) 
head(df.taxa.90)

#####################
####    GAM    ######
#####################

m1 <- gam(Species_richness_5y ~ s(year) + s(Lon, Lat) + s(c.Fy_id) , 
             data = df.taxa.90, family = gaussian, method = 'REML')
m1
summary(m1)
plot(m1)

# with corrected species richness anomaly
m2 <- gam(Species_richness_corr_5y ~ s(year) + s(Lon, Lat) , 
          data = df.taxa.90, family = gaussian, method = 'REML')
m2
summary(m2)
plot(m2)

# model checking
gam.check(m1) # plot fit, convergence with number of iteration, and k check
# check indicated small p-values for s(Lon,Lat): residuals are not randomly distributed
# the residual of the model visualised with histogram looks good
concurvity(m1, full = TRUE) 
concurvity(m1, full = FALSE) 

# specify random effect with gam
# https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/
# make sure to have all random effects as factors
# mgcv can only fit uncorrelated random effects because there’s no way to encode
# a covariance term between the two random effects.
m.re <- gam(Species_richness_5y ~ s(year) + 
              s(study_id, bs="re") + # random intercept
              s(study_id, year, bs="re"), # random slope
              s(c.Fy_id_round, bs="re") + # random intercept
              s(c.Fy_id_round, year, bs="re"), # random slope
         data = df.taxa.90, family = gaussian, method = 'REML')
m.re
summary(m.re)
plot(m.re)

# model checking
gam.check(m.re) # plot fit, convergence with number of iteration, and k check
# check indicated small p-values for s(Lon,Lat): residuals are not randomly distributed
# the residual of the model visualised with histogram looks good
concurvity(m.re, full = TRUE) 
concurvity(m.re, full = FALSE) 

# comparison of models
AIC(m1,m.re) # m1 has lower AIC=128892 compared to m.re=129645

# with corrected species richness anomaly
m.re.corr <- gam(Species_richness_corr_5y ~ s(year) + 
              s(study_id, bs="re") + # random intercept
              s(study_id, year, bs="re"), # random slope
            data = df.taxa.90, family = gaussian, method = 'REML')
m.re.corr
summary(m.re.corr)
plot(m.re.corr)


#####################
####   GAMM4   ######
#####################


# specify random effect with gamm4, see comments by Gavin Simpson here:
# https://stats.stackexchange.com/questions/467234/specifying-random-effects-using-gamm4
# make sure to have all random effects as factors
m.gamm <- gamm4(Species_richness_5y ~ s(year) + 
                  s(study_id, bs="re") + # random intercept
                  s(study_id, year, bs="re"), # random slope
                  s(c.Fy_id_round, bs="re") + # random intercept
                  s(c.Fy_id_round, year, bs="re"), # random slope
            data = df.taxa.90, family = gaussian, REML=TRUE)
m.gamm
summary(m.gamm)
df.taxa.90$predicted1 <- predict(m.gamm$gam)
ggplot(df.taxa.90,aes(x=year,y=predicted1))+geom_point()+facet_wrap(~study_id)+
  geom_line(aes(x=year,y=predicted1),colour="red")
# note: predictions are performed using the observed values of the random variable.
# see Gavin Simpson comments above

# and the same model with lme4 notation (slow run, took 15 min)
# note here I have correlated random intercept and slope
m.gamm4 <- gamm4(Species_richness_5y ~ s(year) ,
                 data = df.taxa.90,
                 random = ~(1+year|c.Fy_id_round) + (1+year|study_id))
m.gamm4
df.taxa.90$predicted2 <- predict(m.gamm4$gam)
ggplot(df.taxa.90,aes(x=year,y=predicted2))+geom_point()+facet_wrap(~study_id)+
  geom_line(aes(x=year,y=predicted2),colour="red")
# note: predictions where the random effects are set to 0.
# see Gavin Simpson comments above



#######################################
#####  GAM Predictive plots     #######
#######################################

# Year partial plot
sort(df.taxa.90$year)
p.year <- seq(min(df.taxa.90$year), max(df.taxa.90$year), by = 1)

newdata = data.frame(year = p.year,
                     Lon = mean(df.taxa.90$Lon), # 10.01777,
                     Lat = mean(df.taxa.90$Lat), # 55.98874, 
                     c.Fy_id = 0
)
head(newdata)
colSums(is.na(newdata)) # no NAs

pgam <- predict(m1, newdata , se.fit = T, type = "response")

# Add the fit and se to the new data frame
newdata$fit <- pgam$fit
newdata$se_lwr <- pgam$fit + qnorm(0.025)* pgam$se.fit # qnorm(0.025) = -1.96
newdata$se_upr <- pgam$fit + qnorm(0.975)* pgam$se.fit # qnorm(0.975) = +1.96


# Plot the fit and se using ggplot2
p.year.plot <- ggplot (newdata, aes (x = year, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "year", y = "richness anomaly") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) +
  ylim(-10, +10) +
  annotate("text", x = 2005, y = 8, label = "all taxonomic ranks", size = 7)

p.year.plot 

ggsave(filename = "~/gam_taxa_richness_anomaly_1990-2020.png", 
       plot = p.year.plot ,
       width = 1539, height = 1300, 
       units = "px"
)



###############################################################
#####  GAM predictive plots with random effect model    #######
###############################################################
library(ggplot2)
# Year partial plot
sort(df.taxa.90$year)
p.year <- seq(min(df.taxa.90$year), max(df.taxa.90$year), by = 1)

newdata = data.frame(year = df.taxa.90$year,
                     study_id = "Sweden_2_1",
                     c.Fy_id_round = "0"
)
head(newdata)
colSums(is.na(newdata)) # no NAs

pgam <- predict(m.re, newdata , se.fit = T, type = "response")

# Add the fit and se to the new data frame
newdata$fit <- pgam$fit
newdata$se_lwr <- pgam$fit + qnorm(0.025)* pgam$se.fit # qnorm(0.025) = -1.96
newdata$se_upr <- pgam$fit + qnorm(0.975)* pgam$se.fit # qnorm(0.975) = +1.96


# Plot the fit and se using ggplot2
p.year.re <- ggplot (newdata, aes (x = year, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "year", y = "richness anomaly") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) +
  ylim(-10, +10)  +
  annotate("text", x = 2000, y = 8, label = "Sweden_2_1", size = 7)

p.year.re 

ggsave(filename = "~/gam_species_richness_anomaly_1990-2020_with_random-effects_Sweden_2_1.png", 
       plot = p.year.re ,
       width = 1539, height = 1300, 
       units = "px"
)


#####################
####    GAMM   ######
#####################
# tutorial https://ge-chunyu.github.io/posts/2024-04-gamm/

# Gavin Simpson:
# assumption of the model is that the observations are conditionally independent
# If you model the autocorrelation through terms in the model, and it is reasonable
# to expect that the smooth functions of Date and the other variables in the model
# are accounting for the temporal structure in the data such that once we consider
# the model, the observations are independent.

# Bart Larsen tutorial
# https://bart-larsen.github.io/GAMM-Tutorial/

# Chunyu Ge tutorial
# https://ge-chunyu.github.io/posts/2024-04-gamm/

# quick look
ggplot(df.taxa.90, aes(x=year, y=Species_richness_corr)) + stat_smooth()

# BEST MODEL STRUCTURE TO LOOK AT TEMPORAL AUTOCORRELATION   
# using the gam notation, working nicely, with uncorrelated random effects
# similar results when s(Lon, Lat) also used as fixed factor, some autocorrelation 
m1 <- gamm(Species_richness_corr_5y ~ s(year) + # s(Lon, Lat) + 
             s(study_id, bs="re") + # random intercept
             s(study_id, year, bs="re") + # random slope
             s(iYear.90, bs="re"), # random intercept to control for non-independence
                               # among samples collected from the same year
          data = df.taxa.90, family = gaussian, method = 'REML')
m1
summary(m1$gam) 

plot(m1$gam)
acf(residuals(m1$gam),main="raw residual ACF", lag.max = 15) 
# autocorrelation exceeding 0.4 at 1 year lag down to 0.2 at 4 year lag

m2 <- gamm(Species_richness_corr_5y ~ s(year) + # s(Lon, Lat) +
             s(study_id, bs="re") + # random intercept
             s(study_id, year, bs="re") + # random slope
             s(iYear.90, bs="re"), # random intercept to control for non-independence
           # among samples collected from the same year
           correlation = corCAR1(form = ~ year|site_id),
           data = df.taxa.90, family = gaussian, method = 'REML')
m2
summary(m2$gam)
plot(m2$gam)
acf(residuals(m2$lme,type="normalized"),main="standardized residual ACF", lag.max = 15)

# Note use of continuous AR1 corCAR1 instead of corAR1 because of the 
# different spacing in year in time series
m2
summary(m2)



