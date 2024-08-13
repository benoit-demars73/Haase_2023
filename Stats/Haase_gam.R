# -------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphics windows 

# -------------------------------------------------------------------
#

library(dplyr) 


setwd("~/Haase_2023/Data")

# retrieve table produce in Taxonomic_richness_anomaly.R script
df.species.90 <- read.table("~/Haase_summary_taxa_rank_species.txt", sep ="\t", header=TRUE)
head(df.species.90)
nrow(df.species.90) # 7495
n_distinct(df.species.90$site_id) # 417
str(df.species.90)
library("lubridate")
df.species.90$site_id <- as.factor(df.species.90$site_id) 

# check for NA values
colSums(is.na(df.species.90))


# Calculate the number of sites for each year
df.site_id_nb <- df.species.90 %>%
  group_by(year) %>%
  summarize(site_id_nb = n_distinct(site_id) )
df.site_id_nb <- as.data.frame(df.site_id_nb)
# calculate sum of species id
df.sum_Fsp_id <- df.species.90 %>%
  group_by(year) %>%
  summarize(sum_Fsp_id = sum(Fsp_id) )
df.sum_Fsp_id <- as.data.frame(df.sum_Fsp_id)
# calculate sum of genus id
df.sum_Fg_id <- df.species.90 %>%
  group_by(year) %>%
  summarize(sum_Fg_id = sum(Fg_id) )
df.sum_Fg_id <- as.data.frame(df.sum_Fg_id)
# calculate sum of family id
df.sum_Ff_id <- df.species.90 %>%
  group_by(year) %>%
  summarize(sum_Ff_id = sum(Ff_id) )
df.sum_Ff_id <- as.data.frame(df.sum_Ff_id)
head(df.sum_Ff_id)

# merge 
df.id <- list(df.site_id_nb, df.sum_Fsp_id, df.sum_Fg_id, df.sum_Ff_id)
df.id <- as.data.frame(df.id)
df.id <- df.id %>%
  select(-c(year.1, year.2, year.3))
head(df.id)
# change to long format
library(tidyr)

df.id.l <- df.id %>% 
  pivot_longer(
    cols = (!c(year, site_id_nb)),
    names_to = "id",
    values_to = "sumF"
  )
df.id.l <- as.data.frame(df.id.l)
head(df.id.l)


library(ggplot2)
# simple bar chart using wide format dataframe
p.sp <- ggplot(data = df.id, aes(x = year, y = site_id_nb)) +
  geom_bar(stat = "identity") +
  labs (x = "year", y = "number of sites") +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) 
p.sp

# stacked bar chart need long format dataframe
p.sp.id <- ggplot(df.id.l, aes(fill=id, y=sumF, x=year)) + 
  geom_bar(position='stack', stat='identity') +
  labs (x = "year", y = "number of sites", fill ="'species'") +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) +
  theme(legend.position = c(.16,.83)) +
  theme(legend.background = element_rect(fill='transparent')) +
  scale_fill_discrete(labels = c("family", "genus", "species")) 
# theme(legend.position = "none") +
# annotate("text", x = 1995, y = 350, label = "'species'", size = 7) 
p.sp.id

ggsave(filename = "~/stacked_bar_chart_species.90_id.png", 
       plot = p.sp.id ,
       width = 1539, height = 1300, 
       units = "px"
)

# stacked bar chart with proportion
p.sp.id.p <- ggplot(df.id.l, aes(fill=id, y=sumF, x=year)) + 
  geom_bar(position='fill', stat='identity') +
  labs (x = "year", y = "Proportion of taxonomic ranks") +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) +
  theme(legend.position = "none") 
p.sp.id.p

ggsave(filename = "~/stacked_bar_chart_species.90_id_proportion.png", 
       plot = p.sp.id.p ,
       width = 1539, height = 1300, 
       units = "px"
)



library(mgcv)
# Noam Ross on spatio-temporal gam https://noamross.github.io/gams-in-r-course/chapter3
# Gavin Simpson on random effect https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/
# bs='re' to have a random effect
summary(df.species.90$Species_richness_5y)
summary(df.species.90$Genus_richness_5y)
summary(df.species.90$Family_richness_5y)

# Trial with richness anomaly
# m.gam.ti <- gam(Species_richness_5y ~ s(year) + s(Lon, Lat) + ti(year,Lon,Lat) + s(c.Fsp_id),
#                data = df.species.90, family = gaussian, method = 'REML')
# m.gam.lm <- gam(Species_richness_5y ~ year + s(Lon, Lat) + s(c.Fsp_id), 
#                data = df.species.90, family = gaussian, method = 'REML') # overall intercept -82 ???

########################
### SPECIES RICHNESS ###

m.gam <- gam(Species_richness_5y ~ s(year) + s(Lon, Lat) + s(c.Fsp_id) , # 
             data = df.species.90, family = gaussian, method = 'REML')


m.gam
summary(m.gam)

# model checking
coef(m.gam) # provide coefficients 
gam.check(m.gam) # plot fit, convergence with number of iteration, and k check
# check indicated small p-values for s(Lon,Lat): residuals are not randomly distributed
# the residual of the model visualised with histogram looks good
concurvity(m.gam, full = TRUE) # overall concurvity (<0.38)
concurvity(m.gam, full = FALSE) # pair-wise concurvity showed ti(year,Lon,Lat) was redundant
AIC(m.gam,m.gam.lm) # model with ti(year,Lon,Lat) had lower AIC, but only + 1% explained deviance
par(mar=c(4,4,1,1))
plot(m.gam,
     se=TRUE, # 95% confidence interval for the mean shape of the effect
     shade=TRUE, # shade.col = "lightblue"
     seWithMean = TRUE, # plot the standard errors of a partial effect term combined with the standard errors of the model intercept
     shift = coef(m.gam)[1], # shift the scale so that the intercept is included
     all.terms=TRUE,
     rug=TRUE,
     residuals=TRUE) # pch = 1, cex = 1
vis.gam(m.gam,
        view = c("Lon", "Lat"),
        plot.type = "contour", 
        contour.col = "black",labcex=0.8,
        too.far = 0.04,
        type = "response",
        main="Richness anomaly spatial structure")

# notes:
# (1) adding country +0.2% explained deviance
# (2) adding study_id +0.5% explained deviance
# (3) so s(Lon, Lat) does its job well !


###################################
#####  Predictive plots     #######
###################################
library(ggplot2)
# Year partial plot
sort(df.species.90$year)
p.year <- seq(min(df.species.90$year), max(df.species.90$year), by = 1)

newdata = data.frame(year = p.year,
                     Lon = mean(df.species.90$Lon), # 10.01777,
                     Lat = mean(df.species.90$Lat), # 55.98874, 
                     c.Fsp_id = mean(df.species.90$c.Fsp_id)
)
head(newdata)
colSums(is.na(newdata)) # no NAs

pgam <- predict(m.gam, newdata , se.fit = T, type = "response")

# Add the fit and se to the new data frame
newdata$fit <- pgam$fit
newdata$se_lwr <- pgam$fit + qnorm(0.025)* pgam$se.fit # qnorm(0.025) = -1.96
newdata$se_upr <- pgam$fit + qnorm(0.975)* pgam$se.fit # qnorm(0.975) = +1.96

head(newdata)
write.table(newdata,"~/species1990-2020_richness-anomaly_fit_417sites.txt",
            sep="\t",row.names=FALSE)

# Plot the fit and se using ggplot2
p.year.plot <- ggplot (newdata, aes (x = year, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "year", y = "richness anomaly") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) +
  ylim(-10, +10) +
  annotate("text", x = 1995, y = 8, label = "'species'", size = 7)

p.year.plot 

ggsave(filename = "~/gam_species1990-2020_richness_anomaly_417sites.png", 
       plot = p.year.plot ,
       width = 1539, height = 1300, 
       units = "px"
)





#########################
###  FAMILY RICHNESS  ###

m.gam <- gam(Family_richness_5y ~ s(year) + s(Lon, Lat) , 
             data = df.species.90, family = gaussian, method = 'REML')


m.gam
summary(m.gam)

# model checking
coef(m.gam) 
gam.check(m.gam) 
concurvity(m.gam, full = TRUE) 
concurvity(m.gam, full = FALSE) 
AIC(m.gam,m.gam.lm) 

#####  Predictive plots     #######

library(ggplot2)
# Year partial plot
sort(df.species.90$year)
p.year <- seq(min(df.species.90$year), max(df.species.90$year), by = 1)

newdata = data.frame(year = p.year,
                     Lon = mean(df.species.90$Lon), # 10.01777,
                     Lat = mean(df.species.90$Lat) # 55.98874, 
                     )
head(newdata)
colSums(is.na(newdata)) # no NAs

pgam <- predict(m.gam, newdata , se.fit = T, type = "response")

# Add the fit and se to the new data frame
newdata$fit <- pgam$fit
newdata$se_lwr <- pgam$fit + qnorm(0.025)* pgam$se.fit # qnorm(0.025) = -1.96
newdata$se_upr <- pgam$fit + qnorm(0.975)* pgam$se.fit # qnorm(0.975) = +1.96

head(newdata)
write.table(newdata,"~/family1990-2020_richness-anomaly_fit_417sites.txt",
            sep="\t",row.names=FALSE)

# Plot the fit and se using ggplot2
p.year.plot <- ggplot (newdata, aes (x = year, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "year", y = "richness anomaly") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) +
  ylim(-10, +10) +
  annotate("text", x = 1995, y = 8, label = "family", size = 7)

p.year.plot 

ggsave(filename = "~/gam_family1990-2020_richness_anomaly_417sites.png", 
       plot = p.year.plot ,
       width = 1539, height = 1300, 
       units = "px"
)



######################################################
########   ALL TAXA 1990-2020   ######################
######################################################

# retrieve table produce in Taxonomic_richness_anomaly.R script
df.taxa.90 <- read.table("~/Haase_summary_taxa_rank_5y_1990-2020.txt", sep ="\t", header=TRUE)
nrow(df.taxa.90) # 18657
n_distinct(df.taxa.90$site_id) # 1120
head(df.taxa.90)
nrow(df.taxa.90) # 18657
str(df.taxa.90)
library("lubridate")
df.taxa.90$site_id <- as.factor(df.taxa.90$site_id) 

# Calculate the number of sites for each year
df.site_id_nb <- df.taxa.90 %>%
  group_by(year) %>%
  summarize(site_id_nb = n_distinct(site_id) )
df.site_id_nb <- as.data.frame(df.site_id_nb)
# calculate sum of species id
df.sum_Fsp_id <- df.taxa.90 %>%
  group_by(year) %>%
  summarize(sum_Fsp_id = sum(Fsp_id) )
df.sum_Fsp_id <- as.data.frame(df.sum_Fsp_id)
# calculate sum of genus id
df.sum_Fg_id <- df.taxa.90 %>%
  group_by(year) %>%
  summarize(sum_Fg_id = sum(Fg_id) )
df.sum_Fg_id <- as.data.frame(df.sum_Fg_id)
# calculate sum of family id
df.sum_Ff_id <- df.taxa.90 %>%
  group_by(year) %>%
  summarize(sum_Ff_id = sum(Ff_id) )
df.sum_Ff_id <- as.data.frame(df.sum_Ff_id)
head(df.sum_Ff_id)

# merge 
df.id <- list(df.site_id_nb, df.sum_Fsp_id, df.sum_Fg_id, df.sum_Ff_id)
df.id <- as.data.frame(df.id)
df.id <- df.id %>%
  select(-c(year.1, year.2, year.3))
head(df.id)
# change to long format
library(tidyr)

df.id.l <- df.id %>% 
  pivot_longer(
    cols = (!c(year, site_id_nb)),
    names_to = "id",
    values_to = "sumF"
  )
df.id.l <- as.data.frame(df.id.l)
head(df.id.l)

library(ggplot2)
# simple bar chart using wide format dataframe
p <- ggplot(data = df.id, aes(x = year, y = site_id_nb)) +
  geom_bar(stat = "identity") +
  labs (x = "year", y = "number of sites") +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) 
  
p
# stacked bar chart need long format dataframe
p.id <- ggplot(df.id.l, aes(fill=id, y=sumF, x=year)) + 
  geom_bar(position='stack', stat='identity') +
  labs (x = "year", y = "number of sites", fill="Rank") +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) +
  theme(legend.position = c(.16,.83)) +
  theme(legend.background = element_rect(fill='transparent')) +
  scale_fill_discrete(labels = c("family", "genus", "species")) 
p.id
ggsave(filename = "~/stacked_bar_chart_taxa_rank_1120sites.png", 
       plot = p.id ,
       width = 1539, height = 1300, 
       units = "px"
)

# stacked bar chart with proportion
p.id.p <- ggplot(df.id.l, aes(fill=id, y=sumF, x=year)) + 
  geom_bar(position='fill', stat='identity') +
  labs (x = "year", y = "Proportion of taxonomic ranks") +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) +
  theme(legend.position = "none") 
p.id.p

ggsave(filename = "~/stacked_bar_chart_taxa_rank_1120sites_proportion.png", 
       plot = p.id.p ,
       width = 1539, height = 1300, 
       units = "px"
)

# gam model
m.gam <- gam(Species_richness_5y ~ s(year) + s(Lon, Lat) + s(c.Fy_id) , 
             data = df.taxa.90, family = gaussian, method = 'REML')


m.gam
summary(m.gam)

# model checking
coef(m.gam) # provide coefficients 
gam.check(m.gam) # plot fit, convergence with number of iteration, and k check
concurvity(m.gam, full = TRUE) # overall concurvity (<0.38)
concurvity(m.gam, full = FALSE) # pair-wise concurvity showed ti(year,Lon,Lat) was redundant
AIC(m.gam,m.gam.lm) 



# save model results
head(df.taxa.90)
df.taxa.90 %>% select(study_id, site_id, Lon, Lat, year, code, Species_richness, Species_richness_5y, c.Fy_id) -> dfx

pp <- predict(m.gam, se.fit = TRUE, type = "response") # since newdata is omitted the predictions are
# based on the data used for the fit.
dfx$fit <- pp$fit
dfx$se_lwr <- pp$fit + qnorm(0.025)* pp$se.fit # qnorm(0.025) = -1.96
dfx$se_upr <- pp$fit + qnorm(0.975)* pp$se.fit # qnorm(0.975) = +1.96

# 95% confidence interval (range or model response uncertainty)
dfx$se_range <- dfx$se_upr - dfx$se_lwr
# glm response residuals
dfx$resp.res <- dfx$Species_richness_5y - dfx$fit
head(dfx)

# write dfx dataframe to file if needed
write.table(dfx,"~/GAM_predictions_1120sites_1990-2020_alltaxa.txt",sep="\t",row.names=FALSE)


# Year partial plot
sort(df.taxa.90$year)
p.year <- seq(min(df.taxa.90$year), max(df.taxa.90$year), by = 1)

newdata = data.frame(year = p.year,
                     Lon = mean(df.taxa.90$Lon), # 10.01777,
                     Lat = mean(df.taxa.90$Lat), # 55.98874, 
                     c.Fy_id = mean(df.taxa.90$c.Fy_id)
)
head(newdata)
colSums(is.na(newdata)) # no NAs

pgam <- predict(m.gam, newdata , se.fit = T, type = "response")

# Add the fit and se to the new data frame
newdata$fit <- pgam$fit
newdata$se_lwr <- pgam$fit + qnorm(0.025)* pgam$se.fit # qnorm(0.025) = -1.96
newdata$se_upr <- pgam$fit + qnorm(0.975)* pgam$se.fit # qnorm(0.975) = +1.96

head(newdata)
write.table(newdata,"~/alltaxa1990-2020_richness_anomaly_fit_1120sites.txt",
            sep="\t",row.names=FALSE)

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

ggsave(filename = "~/gam_species_richness_anomaly_taxa1990-2020.png", 
       plot = p.year.plot ,
       width = 1539, height = 1300, 
       units = "px"
)







########################################################
########   ALL TAXA at FAMILY level id  1990-2020   ####
########################################################

# retrieve table produce in Taxonomic_richness_anomaly.R script
df.family.90 <- read.table("~/Haase_summary_taxa_rank_5y_1990-2020.txt", sep ="\t", header=TRUE)

head(df.family.90)
nrow(df.family.90) # 18657
str(df.family.90)
library("lubridate")
df.family.90$site_id <- as.factor(df.family.90$site_id) 

# gam model
m.gam <- gam(Family_richness_5y ~ s(year) + s(Lon, Lat) , 
             data = df.family.90, family = gaussian, method = 'REML')


m.gam
summary(m.gam)

# model checking
coef(m.gam) # provide coefficients 
gam.check(m.gam) # plot fit, convergence with number of iteration, and k check
concurvity(m.gam, full = TRUE) # overall concurvity (<0.38)
concurvity(m.gam, full = FALSE) # pair-wise concurvity showed ti(year,Lon,Lat) was redundant
AIC(m.gam,m.gam.lm) 



# save model results
head(df.family.90)
df.family.90 %>% select(study_id, site_id, Lon, Lat, year, code, Family_richness, Family_richness_5y) -> dfx

pp <- predict(m.gam, se.fit = TRUE, type = "response") # since newdata is omitted the predictions are
# based on the data used for the fit.
dfx$fit <- pp$fit
dfx$se_lwr <- pp$fit + qnorm(0.025)* pp$se.fit # qnorm(0.025) = -1.96
dfx$se_upr <- pp$fit + qnorm(0.975)* pp$se.fit # qnorm(0.975) = +1.96

# 95% confidence interval (range or model response uncertainty)
dfx$se_range <- dfx$se_upr - dfx$se_lwr
# glm response residuals
dfx$resp.res <- dfx$Family_richness_5y - dfx$fit
head(dfx)

# write dfx dataframe to file if needed
write.table(dfx,"~/GAM_predictions_1120sites_1990-2020_family.txt",sep="\t",row.names=FALSE)


# Year partial plot
sort(df.family.90$year)
p.year <- seq(min(df.family.90$year), max(df.family.90$year), by = 1)

newdata = data.frame(year = p.year,
                     Lon = mean(df.family.90$Lon), # 10.01777,
                     Lat = mean(df.family.90$Lat) # 55.98874, 
                     )
head(newdata)
colSums(is.na(newdata)) # no NAs

pgam <- predict(m.gam, newdata , se.fit = T, type = "response")

# Add the fit and se to the new data frame
newdata$fit <- pgam$fit
newdata$se_lwr <- pgam$fit + qnorm(0.025)* pgam$se.fit # qnorm(0.025) = -1.96
newdata$se_upr <- pgam$fit + qnorm(0.975)* pgam$se.fit # qnorm(0.975) = +1.96

head(newdata)
write.table(newdata,"~/gam fit tables/family1990-2020_richness_anomaly_fit_1120sites.txt",
            sep="\t",row.names=FALSE)

# Plot the fit and se using ggplot2
p.year.plot.f <- ggplot (newdata, aes (x = year, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "year", y = "richness anomaly") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) +
  ylim(-10, +10) +
  annotate("text", x = 1995, y = 8, label = "family", size = 7)

p.year.plot.f

ggsave(filename = "~/gam_family1990-2020_richness_anomaly_1120sites.png", 
       plot = p.year.plot.f ,
       width = 1539, height = 1300, 
       units = "px"
)

