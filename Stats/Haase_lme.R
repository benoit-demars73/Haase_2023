# -------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphics windows 

# -------------------------------------------------------------------
#
library(dplyr) 
library(lme4)

setwd("~/Haase_2023/Data")

################################################################
#####   LME MODELS WITH TAXONOMIC ID LEVEL  - 1816 sites   #####
################################################################


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
summary(df.taxa$c.Fy_id_round)


############### lme site  ####################

m.null <- lmer(Species_richness ~ (1+cYear|site_id),
               data=df.taxa, 
               REML=FALSE,
               control = lmerControl(optimizer = "bobyqa", 
                                     optCtrl=list(maxfun=1e5)) # to avoid convergence issues
)
m.null
summary(m.null)
m <- lmer(Species_richness ~ cYear + (1+cYear|site_id),
          data=df.taxa,
          REML=FALSE,
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl=list(maxfun=1e5)) 
)
m
summary(m) # Intercept = 27.05, se = 0.25 ; Fixed effect cYear = 0.277 , se = 0.013
anova(m.null, m) # Pr(>Chisq) < 2.2e-16 ***
# quick look at model residuals (26664 samples)
m.residuals <- as.data.frame(resid(m))
head(m.residuals)
write.table(m.residuals,"~/lme_1816sites_model-residuals_site.txt",sep="\t",row.names=TRUE,col.names=TRUE)

# Check normality of model residuals with histograms
# checking the normality of the fixed effect
hist((resid(m) - mean(resid(m))) / sd(resid(m)), freq = FALSE)
curve(dnorm, add = TRUE)
# checking the normality of the random effects (random intercept and slope)
hist((ranef(m)$site_id$`(Intercept)` - mean(ranef(m)$site_id$`(Intercept)`)) / sd(ranef(m)$site_id$`(Intercept)`), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m)$site_id$`cYear` - mean(ranef(m)$site_id$`cYear`)) / sd(ranef(m)$site_id$`cYear`), freq = FALSE)
curve(dnorm, add = TRUE)


ranef(m) # coefficients for random effect
# ranef() gives the conditional modes, that is the difference between the (population-level)
# average predicted response for a given set of fixed-effect values (treatment)
# and the response predicted for a particular individual. 
fixef(m) # coefficients for fixed effect
coef(m) # coefficients for fixed + random effects
# is basically just the value of fixef() applicable to each individual plus the value of ranef()
# NOTE: coef cannot be written as a dataframe, but fixef and ranef can


model_coefficients <- as.data.frame(ranef(m))
write.table(model_coefficients,"~/lme model coefficient tables/lme_1816sites.txt",sep="\t",row.names=TRUE,col.names=TRUE)



############### lme site + study ####################

m.null <- lmer(Species_richness ~ (1+cYear|study_id) + (1+cYear|site_id),
               data=df.taxa, 
               REML=FALSE,
               control = lmerControl(optimizer = "bobyqa", 
                                     optCtrl=list(maxfun=1e5)) # to avoid convergence issues
)
m.null
summary(m.null)
m <- lmer(Species_richness ~ cYear + (1+cYear|study_id) + (1+cYear|site_id),
          data=df.taxa,
          REML=FALSE,
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl=list(maxfun=1e5)) 
)
m
summary(m) # Intercept = 28.55, se = 1.59 ; Fixed effect cYear =0.255  , se = 0.085
anova(m.null, m) # Pr(>Chisq) = 0.004674 **
0.255/28.55 # 0.0089 similar to Haase et al 2023
# Check normality of model residuals with histograms
# checking the normality of the fixed effect
hist((resid(m) - mean(resid(m))) / sd(resid(m)), freq = FALSE)
curve(dnorm, add = TRUE)
# checking the normality of the random effects (random intercept and slope)
hist((ranef(m)$site_id$`(Intercept)` - mean(ranef(m)$site_id$`(Intercept)`)) / sd(ranef(m)$site_id$`(Intercept)`), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m)$site_id$`cYear` - mean(ranef(m)$site_id$`cYear`)) / sd(ranef(m)$site_id$`cYear`), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m)$study_id$`(Intercept)` - mean(ranef(m)$study_id$`(Intercept)`)) / sd(ranef(m)$study_id$`(Intercept)`), freq = FALSE)
curve(dnorm, add = TRUE) # lack of normality 
hist((ranef(m)$study_id$`cYear` - mean(ranef(m)$study_id$`cYear`)) / sd(ranef(m)$study_id$`cYear`), freq = FALSE)
curve(dnorm, add = TRUE)

ranef(m) # coefficients for random effect
fixef(m) # coefficients for fixed effect
coef(m) # coefficients for fixed + random effects

model_coefficients <- as.data.frame(ranef(m))
write.table(model_coefficients,"~/lme_1816sites_site-crossed-study.txt",sep="\t",row.names=TRUE,col.names=TRUE)





############### lme site with cFy_id correction ####################

m.null <- lmer(Species_richness ~ (1+cYear|site_id) + (1+cYear|c.Fy_id_round),
               data=df.taxa, 
               REML=FALSE,
               control = lmerControl(optimizer = "bobyqa", 
                                     optCtrl=list(maxfun=1e5)) # to avoid convergence issues
)
m.null
summary(m.null)
m <- lmer(Species_richness ~ cYear + (1+cYear|site_id) + (1+cYear|c.Fy_id_round),
          data=df.taxa,
          REML=FALSE,
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl=list(maxfun=1e5)) 
)
m
summary(m) # Intercept = 26.00, se = 1.95; Fixed effect cYear = 0.283, se = 0.053
anova(m.null, m) # Pr(>Chisq) = 0.0004892 ***

# Check normality of model residuals with histograms
# checking the normality of the fixed effect
hist((resid(m) - mean(resid(m))) / sd(resid(m)), freq = FALSE)
curve(dnorm, add = TRUE)
# checking the normality of the random effects (random intercept and slope)
hist((ranef(m)$site_id$`(Intercept)` - mean(ranef(m)$site_id$`(Intercept)`)) / sd(ranef(m)$site_id$`(Intercept)`), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m)$site_id$`cYear` - mean(ranef(m)$site_id$`cYear`)) / sd(ranef(m)$site_id$`cYear`), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m)$c.Fy_id_round$`(Intercept)` - mean(ranef(m)$c.Fy_id_round$`(Intercept)`)) / sd(ranef(m)$c.Fy_id_round$`(Intercept)`), freq = FALSE)
curve(dnorm, add = TRUE) # lack of normality 
hist((ranef(m)$c.Fy_id_round$`cYear` - mean(ranef(m)$c.Fy_id_round$`cYear`)) / sd(ranef(m)$c.Fy_id_round$`cYear`), freq = FALSE)
curve(dnorm, add = TRUE)

ranef(m) # coefficients for random effect
fixef(m) # coefficients for fixed effect
coef(m) # coefficients for fixed + random effects

model_coefficients <- as.data.frame(ranef(m))
write.table(model_coefficients,"~/lme_1816sites_cFy-site-level.txt",sep="\t",row.names=TRUE,col.names=TRUE)


################################################
#### prepare file for autocorrelation study ####
################################################

# in Excel or Notepad: shift header to the right, insert id column, change grp to site_id header 
df.geo.site <- read.table("site_georeference.txt", sep ="\t", header=TRUE)
head(df.geo.site)
str(df.geo.site$site_id)
df.lme1 <- read.table("~/lme_1816sites_cFy-site-level.txt", sep ="\t", header=TRUE)
df.lme1$site_id <- as.integer(df.lme1$site_id)
df.lme1 <- df.lme1 %>% rename(Trend_est = condval)
head(df.lme1)
nrow(df.lme1) # 3668
df.lme1 <- subset(df.lme1, term %in% c('cYear'))
df.lme1 <- subset(df.lme1, grpvar %in% c('site_id'))
nrow(df.lme1) # 1816
df_merge <- left_join(df.lme1, df.geo.site, by = c("site_id"))  
head(df_merge)
df_merge %>% select(site_id, Trend_est, Lon, Lat) -> df.geo
head(df.geo)
write.table(df.geo,"~/LME_Trend_estimates_1816sites_cFy-site-level_alltaxa.txt",sep="\t",row.names=FALSE)


