#### IN PROGRESS ####

#_____________________________________________
# Exercise 1: Amphibians in SE Wyoming
#_____________________________________________
## Background ##
# Goals of Exercise:
# - Learn R Packages needed for spatial data
# - How to import spatial data in R
# - How to visualize spatial data in R
# - How to analyze spatial data in R

# Data Source: 
# - Amphibian Data: ML Torres
# - Pond Data: ML Torres and USFWS National Wetlands Inventory (website: https://www.fws.gov/program/national-wetlands-inventory)
# - Analyses: ML Torres, MA Murphy, JS Evans, S Albeke

#_____________________________________________
#### 0. Packages ####
#_____________________________________________
install_and_load <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
}

install_and_load(c("sf", 
                   "terra", 
                   "Metrics", 
                   "ROCR", 
                   "spatialEco", 
                   "tidyverse", 
                   "randomForest", 
                   "rfUtilities", 
                   "mapview", 
                   "tmap", 
                   "ggspatial"))
#_____________________________________________
#### 1. Loading Spatial Data ####
#_____________________________________________
##### 1a. Downloading Data #####
# Visit the Google Drive Link below to download the data 
# https://drive.google.com/drive/folders/1esNDgujW7pequkRz7v67DCR4tFlMPhYo?usp=sharing

##### 1b. Setting Working Directory #####
# modify path to the folder you are housing the data. This makes life easier when trying to access the data
path = "C:/Users/melan/OneDrive/Desktop/drive/rworkshop_mlt" 
setwd(path)

##### 1c. Raster Data #####
# let's read in the raster data
# We'll pulling up all the files with the .tif extension
files <- list.files(pattern = "*.tif$", full.names = T)

# creating a raster stack (i.e., a list of all rasters) using `terra::rast`
r <- rast(files)

# plotting the rasters
plot(r)

# Oh no! Our surface relief ratio (SRR) raster's name did not carry over into R correctly.
# Let's name it so we know what this variable is.
names(r[[4]]) <- "srr"

##### 1d. Vector Data #####
# Now, let's read in our vector data
# We have three vectors in the form of CSV files and a shapefile
# As a reminder, shapefiles are spatial files that are *always* vector data (points, lines, polygons)
# First, let's read our CSV files

data <- read.csv("./data.csv")
sites <- read.csv("./SEWY_Sites.csv", header = T)

# the "." = current working directory
# Let's peak at our CSV data:
View(data)
View(sites)

# Note that the `data` file is our amphibian data, whereas the `sites` file is a data frame with spatial data (X and Y) in it.
# However, as a CSV file, `sites` is not a spatial file. There is no datum, and thus no spatial information, associated with this file
# To make this CSV file into a spatial file, *we need to know what datum the data was collected in.* 
# Conveniently, I collected this data - each site's GPS point is in WGS84 format.

sites <- sites %>% 
  st_as_sf(coords = c("X", "Y"), crs = "epsg:4326")
# EPSG = a unique code that all datums have 
# EPSG 4326 is WGS84, a global datum

# Now that sites is a spatial file, let's join the `sites` information to the `data`. 
# Since `data` has a column with the site information in it, we can left join these two data frames by site
df <- data %>% 
  left_join(sites, by = c("SiteName"))

View(df) # take a peak at df and what it looks like!

# To check whether every site was joined correctly, we can search for NAs in the geometry field:
which(is.na(df$geometry), arr.ind = T)

# no NAs exist in this column :)
# If NAs existed, then ensure all site names (or whichever field you're joining your data with) are named correctly!

# Now, let's read in our shapefile. We'll use the `sf` package for this:
wy <- st_read(".", layer = "wy_county")

# let's look at the coordinate reference system (CRS), or datum, of our shp
st_crs(wy)

# Huh. The Wyoming shapefile is in a local datum, NAD83, Zone 13.
# Our spatial data, `df`, is the global datum, WGS84. 
# Not only is NAD83, Zone 13, a local datum (and better captures the geoid in this region), NAD83's horizontal datum is in meters. This datum is better at calculating distance metrics. 
# WGS84, or the typical long/lat we're used to, uses the Earth's gravitational and magnetic model - it's not a great way to measure distances between points. 
# Let's convert our `df` to NAD83, Zone 13:

df <- df %>% 
  st_transform(crs = "epsg:26913")

# or st_transform(crs = st_crs(wy))
# epsg:26913 = NAD83, Zone 13

# Now that we're in the same spatial reference system, let's plot using mapview!
mapview(df, zcol = "Species", col.regions = brewer.pal(6, "Dark2")) +
  mapview(wy, zcol = "COUNTYNAME", col.regions = brewer.pal(23, "Spectral"), legend = F)

#_____________________________________________
#### 2. Analyzing Spatial Data ####
#_____________________________________________
##### 2a. Calculating Landscape Metrics #####
# Landscape Metrics (LM) are the composition and configuration of different land uses within a region
# LMs are important for assessing habitat fragmentation/degradation, connectivity, and more
# We'll manipulate a raster layer, the National Land Cover Database (NLCD), to create separate raster layers important for amphibian presence.
# NLCD is unique as a raster: it contains discrete data classifications for different land uses within each cell (30m)
# While not perfect, NLCD data can help us calculate general trends in land use across North America

# First, let's look at the NLCD layer itself. 
r[[3]]
# classifying the variables:
# 41, 42, 43  forest (frt): mixed, deciduous, evergreen
# 71          open land (old): shrub/scrub
# 81, 82      ag land (ald): hay/pasture, cultivated crops
# 90, 95      wetlands (wet): woody wetlands, emerging herbaceous
# 11          open water (owa): open water

# creating matrices for reclassifying the NLCD data
mfrt <- matrix(c(0, 40, 0,
                 41, 43, 1, 
                 45, 99, 0), ncol = 3, byrow = T)
mold <- matrix(c(0, 70, 0,
                 71, 72, 1,
                 73, 99, 0), ncol = 3, byrow = T)
mald <- matrix(c(0, 80, 0,
                 81, 82, 1,
                 83, 99, 0), ncol = 3, byrow = T)
mwet <- matrix(c(0, 89, 0,
                 90, 95, 1,
                 96, 99, 0), ncol = 3, byrow = T)
mowa <- matrix(c(0, 10, 0,
                 11, 12, 1,
                 13, 99, 0), ncol = 3, byrow = T)

# classifying using terra::classify
frt <- classify(r[[3]], mfrt, right = NA)
old <- classify(r[[3]], mold, right = NA)
ald <- classify(r[[3]], mald, right = NA)
wet <- classify(r[[3]], mwet, right = NA)
owa <- classify(r[[3]], mowa, right = NA)

names(frt)<- "forest"
names(old)<- "open land"
names(ald)<- "ag land"
names(wet)<- "wetlands"
names(owa)<- "open waters"

rlm <- c(frt, old, ald, wet, owa)  

plot(rlm)






