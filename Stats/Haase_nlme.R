# -------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphics windows 

# -------------------------------------------------------------------
#
# general overview:
# https://bbolker.github.io/mixedmodels-misc/notes/corr_braindump.html
# code for nlme:
# https://stackoverflow.com/questions/49796175/how-to-check-and-control-for-autocorrelation-in-a-mixed-effect-model-of-longitud
# Barry McDonald: models with autoregressive errors
# https://r-resources.massey.ac.nz/161251/Lectures/Lecture34.html
# Gavin Simpson on autocorrelation and gam
# https://stats.stackexchange.com/questions/80823/do-autocorrelated-residual-patterns-remain-even-in-models-with-appropriate-corre
# drawing ACF plots with ggplot
# https://stackoverflow.com/questions/72455365/acf-plot-in-ggplot2-from-a-mixed-effects-model
library(dplyr) 
###########################################################
#####   LME MODELS WITH TAXONOMIC ID LEVEL            #####
###########################################################
setwd("~/Haase_2023/Data")
df.taxa <- read.table("Haase_summary_taxa_rank.txt", sep ="\t", header=TRUE)
head(df.taxa)
nrow(df.taxa)
str(df.taxa)
library("lubridate")
df.taxa$sample_id <- as.factor(df.taxa$sample_id)
df.taxa$site_id <- as.factor(df.taxa$site_id) 
df.taxa$country <- as.factor(df.taxa$country) 
df.taxa$study_id <- as.factor(df.taxa$study_id) 
df.taxa$sp_det_prop_round <- as.factor(df.taxa$sp_det_prop_round)
df.taxa$sp_det_factor <- as.factor(df.taxa$sp_det_factor)
df.taxa$c.sp_det_prop_round <- as.factor(df.taxa$c.sp_det_prop_round)
df.taxa$Fsp_id_round <- as.factor(df.taxa$Fsp_id_round)
df.taxa$Fg_id_round <- as.factor(df.taxa$Fg_id_round)
df.taxa$Fm_id_round <- as.factor(df.taxa$Fm_id_round)
df.taxa$Ff_id_round <- as.factor(df.taxa$Ff_id_round)
df.taxa$Fy_id_round <- as.factor(df.taxa$Fy_id_round)
df.taxa$c.Fsp_id_round <- as.factor(df.taxa$c.Fsp_id_round)
df.taxa$c.Fg_id_round <- as.factor(df.taxa$c.Fg_id_round)
df.taxa$c.Fm_id_round <- as.factor(df.taxa$c.Fm_id_round)
df.taxa$c.Ff_id_round <- as.factor(df.taxa$c.Ff_id_round)
df.taxa$c.Fy_id_round <- as.factor(df.taxa$c.Fy_id_round)


n_distinct(df.taxa$site_id) # 1816


# Linear mixed effect model
library(nlme)
m_lme <- lme(Species_richness ~ cYear,
             data=df.taxa,
             random = ~ 1 + cYear | site_id)
m_lme
summary(m_lme)
coef(m_lme)
ranef(m_lme)
# explore autocorrelation

p.lme <- plot(ACF(m_lme),alpha=0.05)
p.lme
ACF(m_lme)
# check normality of model residuals
# with Q-Q plots
qqnorm(resid(m_lme), main = "Q-Q plot for the fixed effect")
qqline(resid(m_lme))
qqnorm(ranef(m_lme)$`(Intercept)`, main = "Q-Q plot for the random intercept")
qqline(ranef(m_lme)$`(Intercept)`)
qqnorm(ranef(m_lme)$`cYear`, main = "Q-Q plot for the random slope")
qqline(ranef(m_lme)$`cYear`)
# histogram better display the density of data points
hist((resid(m_lme) - mean(resid(m_lme))) / sd(resid(m_lme)), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m_lme)$`(Intercept)` - mean(ranef(m_lme)$`(Intercept)`)) / sd(ranef(m_lme)$`(Intercept)`), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m_lme)$`cYear` - mean(ranef(m_lme)$`cYear`)) / sd(ranef(m_lme)$`cYear`), freq = FALSE)
curve(dnorm, add = TRUE)

# if necessary, add temporal correlation to the model 
m_lme_acor <- update(m_lme, correlation = corCAR1(form = ~ iYear|site_id))
m_lme_acor
summary(m_lme_acor)


# to check the ACF after updating the model to include the correlation argument:
# Note: residuals must include the correlation term, so need the normalized residuals,
# the response residuals only relate to the fixed effects part of the model
plot(ACF(m_lme_acor, resType = "normalized"),alpha=0.05)
ACF(m_lme_acor)
# check normality of model residuals 
hist((resid(m_lme_acor) - mean(resid(m_lme_acor))) / sd(resid(m_lme_acor)), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m_lme_acor)$`(Intercept)` - mean(ranef(m_lme_acor)$`(Intercept)`)) / sd(ranef(m_lme_acor)$`(Intercept)`), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m_lme_acor)$`cYear` - mean(ranef(m_lme_acor)$`cYear`)) / sd(ranef(m_lme_acor)$`cYear`), freq = FALSE)
curve(dnorm, add = TRUE)

# compare models with and without AR
AIC(m_lme,m_lme_acor) 
# AIC m_lme       = 183371
# AIC m_lme_acor  = 181810
# model with temporal autocorrelation was statistically better



# The NLME package in R is very powerful for fitting multilevel nonlinear mixed-effects models
# with nested random effects, but it does not fit nonlinear mixed-effects models with 
# crossed random effects (Pinheiro and Bates, 2000)



#####################################################################################



# another idea is just to provide corrected Species_richness as follows (see Fig.4)
df.taxa$Species_richness_corr = df.taxa$Species_richness-(df.taxa$c.Fy_id*12.85)
head(df.taxa)
# Linear mixed effect model
library(nlme)
m_lme <- lme(Species_richness_corr ~ cYear,
             data=df.taxa,
             method="REML",
             random = ~ 1 + cYear | site_id)
m_lme
summary(m_lme)
coef(m_lme)
ranef(m_lme)
# explore autocorrelation
p.lme <- plot(ACF(m_lme),alpha=0.05)
p.lme
# check normality of model residuals 
hist((resid(m_lme) - mean(resid(m_lme))) / sd(resid(m_lme)), freq = FALSE)
curve(dnorm, add = TRUE)
# if necessary, add correlation to the model 
m_lme_acor <- update(m_lme, correlation = corAR1(form = ~ iYear|site_id))
m_lme_acor
summary(m_lme_acor)

# to check the ACF after updating the model to include the correlation argument:
plot(ACF(m_lme_acor, resType = "normalized"),alpha=0.05)
# check normality of model residuals 
hist((resid(m_lme_acor) - mean(resid(m_lme_acor))) / sd(resid(m_lme_acor)), freq = FALSE)
curve(dnorm, add = TRUE)


# compare models with and without AR
AIC(m_lme,m_lme_acor) 
# AIC m_lme       = 181919
# AIC m_lme_acor  = 180766
# model with temporal autocorrelation was statistically better


