# -------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphics windows 

# -------------------------------------------------------------------
#
library(dplyr)
library(nlme)
library (lmtest)
library(ggplot2)

setwd("C:/Users/BDE/OneDrive - NIVA/NIVA/Macroinvertebrate database/Haase paper")
df <- read.table("Manuscript draft/gam fit tables/gam fit table_1990-2020.txt", sep ="\t", header=TRUE)
head(df)

# linear regression 
m <- lm(alltaxa_1120 ~ abundance_1120,data=df)
summary(m)

# Sokal & Rohlf 1995, p.394
# calculate squared difference in residuals 
df.res <- as.data.frame(residuals(m))
df.res <- df.res %>% rename_at('residuals(m)', ~'res')
df.res$year <- seq(1, nrow(df.res))
df.res <- df.res %>% relocate(year, .before = res)
df.res <- df.res %>%
  mutate(dev.sq = (res-mean(df.res$res))^2)
df.res <- df.res %>%
  mutate(diff_res = (lead(res)-res)^2)
head(df.res)
# quick plot
ggplot(df.res, aes(x = df.res$year, y = df.res$diff_res)) +
  geom_point() +
  xlab("year") +
  ylab("diff")
sum(df.res$diff_res, na.rm = TRUE)/sum((df.res$dev.sq)) # 0.34
# < 2 indicate sequences of similar variates, thus autocorrelation

# Durbin Watson autocorrelation test (using year lag = 1), same as Sokal & Rohlf 1995 
# The Durbin Watson statistic ranges from 0 to 4.
# 1.5 < DW < 2.5 indicates no autocorrelation.
# DW < 1.5 indicate positive autocorrelation.
# DW > 2.5 indicate negative autocorrelation.
dwtest(m)
# DW = 0.34, p-value = 1.42e-10, indicating positive serial autocorrelation

# Breusch-Godfrey autocorrelation test (using a range of lags)
bgtest(m, order=3)
# LM test = 22.784, df = 1, p-value = 4.48e-05, indicating significant autocorrelation


# quick model check
plot(m)
# quick look at model residuals time series 
plot(residuals(m)) # some definite trends ....
# autocorrelation plot
acf(residuals(m)) # some strong autocorrelation


# model with AR term
m.ac <- gls(alltaxa_1120 ~ abundance_1120, data=df, 
              correlation = corAR1(form=~year),
              na.action=na.omit)
summary(m.ac)
# autocorrelation plot
plot(ACF(m.ac, resType = "normalized"),alpha=0.05)
# still lots of autocorrelation ... so AR1 does not seem sufficient ...





# Let's try to thin the data by removing odd years
df.thin <- df %>% filter(year %% 2 == 0)
head(df.thin)

m.thin <- lm(alltaxa_1120 ~ abundance_1120,data=df.thin)
summary(m.thin)
# quick model check
plot(m.thin)
# quick look at model residuals time series 
plot(residuals(m.thin)) # some apparent trends ....
# autocorrelation plot
acf(residuals(m.thin)) # autocorrelation reduced to within 95% confidence interval

# Durbin Watson autocorrelation test (using year lag = 1)
dwtest(m.thin)
# DW = 1.0168, p-value = 0.007352, indicating positive serial autocorrelation

# Breusch-Godfrey autocorrelation test (using a range of lags)
bgtest(m.thin, order=3)
# LM test = 7.47, df = 1, p-value = 0.05842, indicating no significant autocorrelation



# removing year 1990, 1991, 2020 may help 
df.92 <- df %>% filter(!year %in% c("1990", "1991", "2020"))
head(df.92)

m.92 <- lm(alltaxa_1120 ~ abundance_1120,data=df.92)
summary(m.92)
# quick model check
plot(m.92)
# quick look at model residuals time series 
plot(residuals(m.92)) # some apparent trends ....
# autocorrelation plot
acf(residuals(m.92)) # significant autocorrelation 

# Durbin Watson autocorrelation test
dwtest(m.92)
# DW = 0.36787, p-value = 1.822e-09, indicating positive serial autocorrelation

# Breusch-Godfrey autocorrelation test
bgtest(m.92, order=3)
# LM test = 23.8, df = 1, p-value = 2.8e-05, indicating significant autocorrelation

# need thinning as well, so let's just use full dataset thinned by eliminating odd years



###########################
###  family 1120 sites  ###
###########################

m.thin <- lm(family_1120 ~ abundance_1120,data=df.thin)
summary(m.thin)
# quick model check
plot(m.thin)
# quick look at model residuals time series 
plot(residuals(m.thin)) # some apparent trends ....
# autocorrelation plot
acf(residuals(m.thin)) # autocorrelation reduced to within 95% confidence interval

# Durbin Watson autocorrelation test (using year lag = 1)
dwtest(m.thin)
# DW = 0.79, p-value = 0.0011, indicating positive serial autocorrelation

# Breusch-Godfrey autocorrelation test (using a range of lags)
bgtest(m.thin, order=3)
# LM test = 5.31, df = 1, p-value = 0.15, indicating no significant autocorrelation



#############################
###  'species' 417 sites  ###
#############################

df.thin <- df %>% filter(year %% 2 == 0)
head(df.thin)


m.thin <- lm(alltaxa_417 ~ abundance_417,data=df.thin)
summary(m.thin)
# quick model check
plot(m.thin)
# quick look at model residuals time series 
plot(residuals(m.thin)) # some apparent trends ....
# autocorrelation plot
acf(residuals(m.thin)) # autocorrelation reduced to within 95% confidence interval

# Durbin Watson autocorrelation test (using year lag = 1)
dwtest(m.thin)
# DW = 0.57, p-value = 0.00018, indicating positive serial autocorrelation

# Breusch-Godfrey autocorrelation test (using a range of lags)
bgtest(m.thin, order=3)
# LM test = 5.31, df = 1, p-value = 0.1006, indicating significant autocorrelation


# model with AR term
m.ac <- gls(alltaxa_417 ~ abundance_417, data=df.thin, 
            correlation = corAR1(form=~year),
            na.action=na.omit)
summary(m.ac)
# autocorrelation plot
plot(ACF(m.ac, resType = "normalized"),alpha=0.05)


#############################
###   family 417 sites    ###
#############################

df.thin <- df %>% filter(year %% 2 == 0)
head(df.thin)


m.thin <- lm(family_417 ~ abundance_417,data=df.thin)
summary(m.thin)
# quick model check
plot(m.thin)
# quick look at model residuals time series 
plot(residuals(m.thin)) # some apparent trends ....
# autocorrelation plot
acf(residuals(m.thin)) # autocorrelation reduced to within 95% confidence interval

# Durbin Watson autocorrelation test (using year lag = 1)
dwtest(m.thin)
# DW = 0.21, p-value = 3.412e-08, indicating positive serial autocorrelation

# Breusch-Godfrey autocorrelation test (using a range of lags)
bgtest(m.thin, order=3)
# LM test = 10.301, df = 3, p-value = 0.016, indicating significant autocorrelation


# model with AR term
m.ac <- gls(family_417 ~ abundance_417, data=df.thin, 
            correlation = corAR1(form=~year),
            na.action=na.omit)
summary(m.ac)
# autocorrelation plot
plot(ACF(m.ac, resType = "normalized"),alpha=0.05)
