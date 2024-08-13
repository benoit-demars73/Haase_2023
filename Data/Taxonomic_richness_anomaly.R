# -------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphics windows 

# -------------------------------------------------------------------
#
library(dplyr) 
library(readr) # read.delim function

setwd("~/Haase_2023/Data")


######################################################
####  subset time series with richness ANOMALY    ####
######################################################

my_data.xy <- read.table("Haase_summary_taxa_rank.txt", sep ="\t", header=TRUE)
head(my_data.xy)
# Calculate the five years (2003-2007) reference average richness per site 
# anomaly calculations
# https://stackoverflow.com/questions/16420267/calculation-of-anomalies-on-time-series


# SPECIES richness
my_data.xy_5y <- transform(my_data.xy, 
                           Species_richness_5y=Species_richness - ave(replace(Species_richness, year < 2003 | year > 2007, NA), 
                                                                      site_id, FUN=function(t) mean(t, na.rm=TRUE)))
head(my_data.xy_5y)
sum(is.na(my_data.xy_5y$Species_richness_5y)) # 2158, about 8% data loss
2158/26664

# remove records with no richness average during 2003-2007 period
my_data.xy_5y <- my_data.xy_5y %>% filter(!is.na(Species_richness_5y))
n_distinct(my_data.xy_5y$site_id) # 1598

# Calculate the range of time series 
site_range_year <- my_data.xy_5y %>%
  group_by(site_id) %>%
  summarize(range_year = max(year)-min(year)) 
# Calculate the number of time points in the time series 
site_n_year <- my_data.xy_5y %>%
  group_by(site_id) %>%
  summarize(count_year = n_distinct(year)) 
# check out how many sites were removed in previous step (5 years reference period)
nrow(site_range_year) # so 12% of the 1816 sites in Haase et al (2023)
1598/1816 # 88% 

# Merge the year range back into the original dataframe
df.range <- left_join(my_data.xy_5y, site_range_year, by = "site_id")
head(df.range)
# Merge the number of time points into the original dataframe
df.range.point <- left_join(df.range, site_n_year, by = "site_id")
head(df.range.point)
# subset site with minimum x time points and min x year range
df.subset1 <- df.range.point[df.range.point$range_year >= 15 & 
                               df.range.point$count_year >= 8, ]

nrow(my_data.xy) # 26664
nrow(df.range.point) # 24506
nrow(df.subset1) # 19911
19911/26664 # 75%
n_distinct(my_data.xy$site_id) # 1816
n_distinct(df.range.point$site_id) # 1598
n_distinct(df.subset1$site_id) # 1191
1191/1816 # 66%



head(df.subset1)


# GENUS richness
my_data_5y_gen <- transform(df.subset1, 
                            Genus_richness_5y=Genus_richness - ave(replace(Genus_richness, year < 2003 | year > 2007, NA), 
                                                                   site_id, FUN=function(t) mean(t, na.rm=TRUE)))
head(my_data_5y_gen)
sum(is.na(my_data_5y_gen$Genus_richness_5y)) # 2158, about 8% data loss
2158/26668

# remove records with no richness average during 2003-2007 period
my_data_5y_gen <- my_data_5y_gen %>% filter(!is.na(Genus_richness_5y))
n_distinct(my_data_5y_gen$site_id) # 1191




# FAMILY richness
my_data_5y_fam <- transform(my_data_5y_gen, 
                            Family_richness_5y=Family_richness - ave(replace(Family_richness, year < 2003 | year > 2007, NA), 
                                                                     site_id, FUN=function(t) mean(t, na.rm=TRUE)))
head(my_data_5y_fam)
sum(is.na(my_data_5y_fam$Family_richness_5y)) # 2158, about 8% data loss
2158/26668

# remove records with no richness average during 2003-2007 period
my_data_5y_fam <- my_data_5y_fam %>% filter(!is.na(Family_richness_5y))
n_distinct(my_data_5y_fam$site_id) # 1191




# write file
write.table(my_data_5y_fam,"Haase_summary_taxa_rank_5y.txt",
            sep="\t",row.names=FALSE)




################################################################
####  subset time series with richness anomaly for species  ####
################################################################

my_data.xy <- read.table("Haase_summary_taxa_rank_5y.txt", sep ="\t", header=TRUE)
head(my_data.xy)
str(my_data.xy)
library("lubridate")
my_data.xy$site_id <- as.factor(my_data.xy$site_id) 

str(my_data.xy$site_id) 
nrow(my_data.xy) # 19911

# subset samples with identification to species
df.species <- subset(my_data.xy, TaxonomicRes %in% c('species'))
nrow(df.species) # 8624
8624/11861 # 73% of original data records
n_distinct(df.species$site_id) # 481

# sample subset 1990-2019
df.species.90 <- df.species[df.species$year > 1989 & df.species$year < 2021, ]
nrow(df.species.90) # 7963
# Calculate the range of time series 
site_range_year <- df.species.90 %>%
  group_by(site_id) %>%
  summarize(range_year_sp = max(year)-min(year)) 
# Calculate the number of time points in the time series 
site_n_year <- df.species.90 %>%
  group_by(site_id) %>%
  summarize(count_year_sp = n_distinct(year)) 
# Merge the year range back into the original dataframe
df.range <- left_join(df.species.90, site_range_year, by = "site_id")
head(df.range)
# Merge the number of time points into the original dataframe
df.range.point <- left_join(df.range, site_n_year, by = "site_id")
head(df.range.point)
# subset site with minimum 8 time points and min 15 year range
df.subset1 <- df.range.point[df.range.point$range_year_sp >= 15 & 
                               df.range.point$count_year_sp >= 8, ]

nrow(my_data.xy) # 19911
nrow(df.range.point) # 7963
nrow(df.subset1) # 7495
7495/19911 # 37%
n_distinct(my_data.xy$site_id) # 1191
n_distinct(df.range.point$site_id) # 481
n_distinct(df.subset1$site_id) # 417
417/481 # 87%

head(df.subset1)
# check for NA values
colSums(is.na(df.subset1))

# write file
write.table(df.subset1,"Haase_summary_taxa_rank_species.txt",
            sep="\t",row.names=FALSE)





############################################################################
####  1990-2020 subset time series with richness anomaly for all taxa   ####
############################################################################


my_data.xy <- read.table("Haase_summary_taxa_rank_5y.txt", sep ="\t", header=TRUE)
head(my_data.xy)
str(my_data.xy)
library("lubridate")
my_data.xy$site_id <- as.factor(my_data.xy$site_id) 

str(my_data.xy$site_id) 
nrow(my_data.xy) # 19911

# sample subset 1990-2019
df.taxa.90 <- my_data.xy[my_data.xy$year > 1989 & my_data.xy$year < 2021, ]
nrow(df.taxa.90) # 19178
# Extract first year of time series 
site_first_year <- df.taxa.90 %>%
  group_by(site_id) %>%
  summarize(first_year = min(year)) 
# Extract last year of time series 
site_last_year <- df.taxa.90 %>%
  group_by(site_id) %>%
  summarize(last_year = max(year)) 
# Calculate the range of time series 
site_range_year <- df.taxa.90 %>%
  group_by(site_id) %>%
  summarize(range_year_sp = max(year)-min(year)) 
# Calculate the number of time points in the time series 
site_n_year <- df.taxa.90 %>%
  group_by(site_id) %>%
  summarize(count_year_sp = n_distinct(year)) 
# Merge the year range back into the original dataframe
df.year1 <- left_join(df.taxa.90, site_first_year, by = "site_id")
head(df.year1)
# Merge the year range back into the original dataframe
df.last.year <- left_join(df.year1, site_last_year, by = "site_id")
head(df.last.year)
# Merge the year range back into the original dataframe
df.range <- left_join(df.last.year, site_range_year, by = "site_id")
head(df.range)
# Merge the number of time points into the original dataframe
df.range.point <- left_join(df.range, site_n_year, by = "site_id")
head(df.range.point)
# subset site with minimum 8 time points and min 15 year range
df.subset2 <- df.range.point[df.range.point$range_year_sp >= 15 & 
                               df.range.point$count_year_sp >= 8, ]

nrow(my_data.xy) # 19911
nrow(df.range.point) # 19178
nrow(df.subset2) # 18657
18657/19911 # 94%
n_distinct(my_data.xy$site_id) # 1191
n_distinct(df.range.point$site_id) # 1191
n_distinct(df.subset2$site_id) # 1120
1120/1191 # 94%

head(df.subset2)
# check for NA values
colSums(is.na(df.subset2))

# write file
write.table(df.subset2,"Haase_summary_taxa_rank_5y_1990-2020.txt",
            sep="\t",row.names=FALSE)



