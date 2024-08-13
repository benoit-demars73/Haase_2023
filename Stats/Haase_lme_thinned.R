# -------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphics windows 

# -------------------------------------------------------------------
#
library(dplyr) 
library("lubridate")
library(lme4)
library(spThin)

setwd("~/Haase_2023/Data")

###########################################################
#####   LME MODELS WITH TAXONOMIC ID LEVEL            #####
###########################################################

df.taxa <- read.table("Haase_summary_taxa_rank.txt", sep ="\t", header=TRUE)
head(df.taxa)
nrow(df.taxa)
str(df.taxa)

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


# subset of sites to avoid spatial autocorrelation

df.geo <- read.table("site_georeference.txt", sep ="\t", header=TRUE)
df.geo$species <- 'sp'
head(df.geo)

thinned <- thin(
  df.geo,
  lat.col = "Lat",
  long.col = "Lon",
  spec.col = "species",
  thin.par = 30,
  reps = 100,
  locs.thinned.list.return = FALSE,
  write.files = TRUE,
  max.files = 5,
  out.dir = "Outputs",
  out.base = "thinned_data",
  write.log.file = TRUE,
  log.file = "spatial_thin_log.txt",
  verbose = TRUE
)

# subset df.taxa
library(tidyverse)
df.thin20 <- read.table("~/thinned_data_20km.txt", sep ="\t", header=TRUE)
head(df.thin20)
df.taxa.20 <- left_join(df.thin20, df.taxa, by = c("Lat" = "Lat", "Lon" = "Lon"))
head(df.taxa.20)

# run lme model
############### lme site  ####################

m.null <- lmer(Species_richness ~ (1+cYear|site_id),
               data=df.taxa.20, 
               REML=FALSE,
               control = lmerControl(optimizer = "bobyqa", 
                                     optCtrl=list(maxfun=1e5)) # to avoid convergence issues
)
m.null
summary(m.null)
m <- lmer(Species_richness ~ cYear + (1+cYear|site_id),
          data=df.taxa.20,
          REML=FALSE,
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl=list(maxfun=1e5)) 
)
m
summary(m) # Intercept = 28.72, se = 0.41 ; Fixed effect cYear = 0.278 , se = 0.022
anova(m.null, m) # Pr(>Chisq) < 2.2e-16 ***

# Check normality of model residuals with histograms
# checking the normality of the fixed effect
hist((resid(m) - mean(resid(m))) / sd(resid(m)), freq = FALSE)
curve(dnorm, add = TRUE)
# checking the normality of the random effects (random intercept and slope)
hist((ranef(m)$site_id$`(Intercept)` - mean(ranef(m)$site_id$`(Intercept)`)) / sd(ranef(m)$site_id$`(Intercept)`), freq = FALSE)
curve(dnorm, add = TRUE)
hist((ranef(m)$site_id$`cYear` - mean(ranef(m)$site_id$`cYear`)) / sd(ranef(m)$site_id$`cYear`), freq = FALSE)
curve(dnorm, add = TRUE)



############### lme site + study ####################

m.null <- lmer(Species_richness ~ (1+cYear|study_id) + (1+cYear|site_id),
               data=df.taxa.20, 
               REML=FALSE,
               control = lmerControl(optimizer = "bobyqa", 
                                     optCtrl=list(maxfun=1e5)) # to avoid convergence issues
)
m.null
summary(m.null)
m <- lmer(Species_richness ~ cYear + (1+cYear|study_id) + (1+cYear|site_id),
          data=df.taxa.20,
          REML=FALSE,
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl=list(maxfun=1e5)) 
)
m
summary(m) # Intercept = 29.33, se = 1.69 ; Fixed effect cYear =0.224  , se = 0.076
anova(m.null, m) # Pr(>Chisq) = 0.008995 **

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



############### lme site with cFy_id correction ####################

m.null <- lmer(Species_richness ~ (1+cYear|site_id) + (1+cYear|c.Fy_id_round),
               data=df.taxa.20, 
               REML=FALSE,
               control = lmerControl(optimizer = "bobyqa", 
                                     optCtrl=list(maxfun=1e5)) # to avoid convergence issues
)
m.null
summary(m.null)
m <- lmer(Species_richness ~ cYear + (1+cYear|site_id) + (1+cYear|c.Fy_id_round),
          data=df.taxa.20,
          REML=FALSE,
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl=list(maxfun=1e5)) 
)
m
summary(m) # Intercept = 24.31, se = 1.59; Fixed effect cYear = 0.306, se = 0.026
anova(m.null, m) # Pr(>Chisq) = 0.0003146 ***

ranef(m) # coefficients for random effect
fixef(m) # coefficients for fixed effect
coef(m) # coefficients for fixed + random effects

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

model_coefficients <- as.data.frame(ranef(m))
write.table(model_coefficients,"~/lme_1816sites_cFy-site-level_thin20km.txt",sep="\t",row.names=TRUE,col.names=TRUE)




##########  lme site + study with c.Fy_id correction  ################

m.null <- lmer(Species_richness ~ (1+cYear|study_id) + (1+cYear|site_id) + 
                 (1+cYear|c.Fy_id_round) ,
               data=df.taxa.20, 
               REML=FALSE,
               control = lmerControl(optimizer = "bobyqa", 
                                     optCtrl=list(maxfun=1e5)) # to avoid convergence issues
)
m.null
summary(m.null)
m <- lmer(Species_richness ~ cYear + (1+cYear|study_id) + (1+cYear|site_id) + 
            (1+cYear|c.Fy_id_round) , 
          data=df.taxa.20,
          REML=FALSE,
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl=list(maxfun=1e5)) 
)
m
summary(m) # Intercept = 24.64, se = 2.25 ; Fixed effect cYear = 0.289 , se = 0.075
anova(m.null, m) # Pr(>Chisq) = 0.003354 **

ranef(m) # coefficients for random effect
fixef(m) # coefficients for fixed effect
coef(m) # coefficients for fixed + random effects


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
hist((ranef(m)$c.Fy_id_round$`(Intercept)` - mean(ranef(m)$c.Fy_id_round$`(Intercept)`)) / sd(ranef(m)$c.Fy_id_round$`(Intercept)`), freq = FALSE)
curve(dnorm, add = TRUE) # lack of normality 
hist((ranef(m)$c.Fy_id_round$`cYear` - mean(ranef(m)$c.Fy_id_round$`cYear`)) / sd(ranef(m)$c.Fy_id_round$`cYear`), freq = FALSE)
curve(dnorm, add = TRUE)


model_coefficients <- as.data.frame(ranef(m))
write.table(model_coefficients,"~/lme_1816sites_cFy_study_level_thin20km.txt",sep="\t",row.names=TRUE,col.names=TRUE)



#######################################################
#####   LME MODELS WITH FAMILY ID LEVEL           #####
#######################################################



############### lme site  ####################

m.null <- lmer(Family_richness ~ (1+cYear|site_id),
               data=df.taxa.20, 
               REML=FALSE,
               control = lmerControl(optimizer = "bobyqa", 
                                     optCtrl=list(maxfun=1e5)) # to avoid convergence issues
)
m.null
summary(m.null)
m <- lmer(Family_richness ~ cYear + (1+cYear|site_id),
          data=df.taxa.20,
          REML=FALSE,
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl=list(maxfun=1e5)) 
)
m
summary(m) # Intercept = 23.20, se= 0.30; Fixed effect cYear = 0.200 , se = 0.016
anova(m.null, m) # Pr(>Chisq) < 2.2e-16 ***

ranef(m) # coefficients for random effect
fixef(m) # coefficients for fixed effect
coef(m) # coefficients for fixed + random effects

model_coefficients <- as.data.frame(ranef(m))
write.table(model_coefficients,"~/lme_Family_1816sites_thin20km.txt",sep="\t",row.names=TRUE,col.names=TRUE)




############### lme site + study  ####################

m.null <- lmer(Family_richness ~ (1+cYear|study_id) + (1+cYear|site_id),
               data=df.taxa.20, 
               REML=FALSE,
               control = lmerControl(optimizer = "bobyqa", 
                                     optCtrl=list(maxfun=1e5)) # to avoid convergence issues
)
m.null
summary(m.null)
m <- lmer(Family_richness ~ cYear + (1+cYear|study_id) + (1+cYear|site_id),
          data=df.taxa.20,
          REML=FALSE,
          control = lmerControl(optimizer = "bobyqa", 
                                optCtrl=list(maxfun=1e5)) 
)
m
summary(m) # Intercept = 21.31, se = 0.97 ; Fixed effect cYear = 0.146, se = 0.034
anova(m.null, m) # Pr(>Chisq) = 0.0006914 ***

ranef(m) # coefficients for random effect
fixef(m) # coefficients for fixed effect
coef(m) # coefficients for fixed + random effects

model_coefficients <- as.data.frame(ranef(m))
write.table(model_coefficients,"~/lme_Family_1816sites_site-crossed-study_thin20km.txt",sep="\t",row.names=TRUE,col.names=TRUE)




