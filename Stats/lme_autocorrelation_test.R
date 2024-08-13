# -------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphics windows 

# -------------------------------------------------------------------
#


library(dplyr) 
library(ncf)

setwd("~/Haase_2023/Data")

###########################################################
########   ALL TAXA 1816 sites with correction   ##########
###########################################################

# retrieve table produced in Haase_lme.R
df.geo <- read.table("~/LME_Trend_estimates_1816sites_cFy-site-level_alltaxa.txt", sep ="\t", header=TRUE)
head(df.geo)

library(sf) 
#--- check the class ---#
class(df.geo) # data.frame
#--- recognize it as an sf ---#
df.geo_sf <- st_as_sf(df.geo, coords = c("Lon","Lat"))
#--- take a look at the data ---#
head(df.geo_sf)
#--- set crs ---#
df.geo_sf <- st_set_crs(df.geo_sf, 4326) # WGS84
#--- see the change ---#
head(df.geo_sf) # 
nrow(df.geo_sf)
# transform coordinate reference system, reprojection
df.geo_3035 <- st_transform(df.geo_sf, crs = 3035) # EU EPSG:3035
head (df.geo_3035) 

# write new geometry in column
df.geo_3035 <- df.geo_3035 %>%   
  mutate(     y = st_coordinates(df.geo_3035)[, 2],
              x = st_coordinates(df.geo_3035)[, 1]   )
head (df.geo_3035)
nrow (df.geo_3035) # 1816
# converting back to a dataframe by dropping the geometry
df.geo_3035_df <- df.geo_3035 %>% st_drop_geometry()
class(df.geo_3035_df)
head (df.geo_3035_df)
nrow(df.geo_3035_df)

df.geo_3035_df <- df.geo_3035_df %>% rename(X_3035 = x,
                                            Y_3035 = y)
head (df.geo_3035_df)

# Transform coordinate units to km
df.geo_3035_df$X_km <- df.geo_3035_df$X_3035/1000 # longitude in km
df.geo_3035_df$Y_km <- df.geo_3035_df$Y_3035/1000 # latitude in km

library(ncf)
richness.corr <- spline.correlog(
  df.geo_3035_df$X_km, # units are km
  df.geo_3035_df$Y_km, # units are km
  df.geo_3035_df$Trend_est,
  w = NULL,
  df = NULL,
  type = "boot",
  resamp = 99,
  npoints = 50, # every 10 km with xmax = 500 km
  save = FALSE,
  filter = FALSE,
  fw = 0,
  max.it = 25,
  xmax = 500, # units are km
  latlon = FALSE,
  na.rm = FALSE,
  quiet = FALSE
)
richness.corr
summary(richness.corr)

autocorr.plot <- plot(richness.corr) # save as .tiff

print.default(richness.corr)

