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
                   "ggspatial", 
                   "RColorBrewer",
                   "landscapemetrics"))
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

# Oh no! We have two unknown rasters in our list! Let's see what they are by looking at the files object:
files

# lyr.1 = hli (heat load index)
# lyr1 = srr (surface relief ratio)
# Also, NED30m_WY (our digital elevation model, DEM) and NLCD Land Cover Class are a mouth full. 
# let's rename these rasters!
nms <- list.files(pattern = "*.tif$") 
nms <- gsub("(\\.tif).*", "", nms)

names(r) <- nms

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

sites <- st_as_sf(sites, coords = c("X", "Y"), crs = "epsg:4326")
# EPSG = a unique code that all datums have 
# EPSG 4326 is WGS84, a global datum

# Let's plot this using mapview!
mapview(sites)

# Now that sites is a spatial file, let's join the `sites` information to the `data`. 
# Since `data` has a column with the site information in it, we can left join these two data frames by site
df <- data %>% 
  left_join(sites, by = c("SiteName")) %>% 
  st_as_sf()

View(df) # take a peak at df and what it looks like!

# To check whether every site was joined correctly, we can search for NAs in the geometry field:
which(is.na(df$geometry), arr.ind = T)

# no NAs exist in this column :)
# If NAs existed, then ensure all site names (or whichever field you're joining your data with) are named correctly!

# Now, let's read in our shapefile. We'll use the `sf` package for this:
wy <- st_read(".", layer = "wy_county")

# let's look at the coordinate reference system (CRS) of wyoming
st_crs(wy)

# note: the WY shapefile's CRS is a local datum, NAD83
# Also note: this file uses the local map projection, Universal Transverse Mercator (UTM), Zone 13
# Our spatial data, df, is the global datum, WGS84. 
# Not only is NAD83 a local datum for North America (and better captures the geoid in this region), and UTM Zone 13 is a projection that is specific for our region in Wyoming, the NAD83 coordinate system is in meters. This datum is better at calculating distance metrics. 
# WGS84, or the typical long/lat we're used to, is a global datum and uses a global projection (WGS 84) that uses the Earth's gravitational and magnetic model - it's not a great way to measure distances between points. 
# let's transform our data frame to NAD83, Zone 13

df <- df %>% 
  st_transform(crs = "epsg:26913")

# ...or st_transform(crs = st_crs(wy))
# epsg:26913 = NAD83, Zone 13

# Now that we're in the same spatial reference system, let's plot using mapview!
mapview(df, zcol = "Species", col.regions = brewer.pal(6, "Dark2")) +
  mapview(wy, zcol = "COUNTYNAME", col.regions = brewer.pal(23, "Spectral"), legend = F)

# So, WY is the entire state of Wyoming, whereas we are interested in shrinking our study area to capture just the points. Let's create a study area shapefile 

# Let's see what the min and max is for our sites:
st_bbox(df)

# let's use the min/max of our points to create a polygon:
study <- st_as_sfc(st_bbox(df))

plot(study)
plot(df$geometry, add = T)

# Note that the bounding box is where the min/max points are. So, we'll increase study area by adding a 1000m buffer around the bounding box to include these points! 
# this also helps with preventing "edge effects" when we work with raster data - we don't want the boundary of our study area to be the reason why our analyses gets messed up

buff <- st_buffer(study, dist = 1000)

# cool! We have a study area! Let's crop our earlier rasters to our buffered study area:
rc <- terra::crop(r, buff)

# Now that our rasters are cropped to size, let's calculate a few more variables, specifically from our elevation layer (dem)
# slope (terra)
slope <- terrain(r["dem"], v = "slope", neighbors = 8, unit = "degrees")

# aspect (terra)
asp <- terrain(r["dem"], v = "aspect", neighbors = 8, unit = "degrees")

# we can also calculate HLI and SRR using the spatialEco package! However, for time's sake, I included those rasters from the beginning :)

# creating a stack of the rasters that have same extent, origin, crs, resolution
rs <- c(r["dem"], r["ffp"], r["gsp"], r["srr"], r["hli"], slope, asp)

# we'll use the nlcd layer from the `r` raster stack in step 2c.

#_____________________________________________
#### 2. Analyzing Spatial Data ####
#_____________________________________________
##### 2a. Creating Pseudoabsences #####
#_____________________________________________
# While I have data for many species, if I analyze it as is, my analyses assumes (a) perfect detection and (b) I located all possible locations for amphibians in the study area. 
# I am good at finding amphibians, but I'm not that good!
# So, we will prepare our data with that in mind!
# let's create pseudoabsences! 

# to start, let's do one species for time's sake - boreal chorus frog (Pseudacris maculata, PSMA)
psma <- df %>% 
  filter(Species == "PSMA")

names(psma)

# let's create a random sample of points for PSMA within the study area
rp <- st_sample(study, size = 2*nrow(psma), type = "random") %>% 
  st_as_sf(data.frame(PointID = 1:nrow(psma), Present = 0)) %>% 
  rename(geometry = x)

# let's add the random points to the presence points
dat.psma <- psma %>% 
  mutate(Present = 1) %>% 
  rename(PointID = LabID) %>% 
  select(PointID, Present, geometry) %>% 
  rbind(rp)

# sometimes, R gets confused, so we're ensuring our dat.psma is still a spatial file
st_geometry(dat.psma) <- "geometry"

# le'ts visualize the data using ggplot!
ggplot() +
  geom_sf(data = buff) +
  geom_sf(data = dat.psma, aes(group=as.factor(Present), col = as.factor(Present)))

#_____________________________________________
#### 2b. Extracting Raster Variables ####
#_____________________________________________
# For modeling purposes, we now need our predictor (independent) variables for PSMA presence.
# So, let's extract data from our rasters!
e <- data.frame(terra::extract(rs, vect(dat.psma)))

# let's look at the results! Mainly, we're looking for any NAs or any outliers in the data to further investigate
summary(e)

which(is.na(dat.psma), arr.ind = T)

# combine covariates with our response variable (presence/absence)

dat.psma <- dat.psma %>% 
  bind_cols(e) %>% 
  mutate(Present = as.factor(Present)) %>% 
  select(-ID)


#_____________________________________________
##### 2c. Calculating Landscape Metrics #####
#_____________________________________________
# Other predictor variables can be calculated using Landscape Metrics
# Landscape Metrics describe the composition and configuration of different land uses within a region
# LMs are important for assessing habitat fragmentation/degradation, connectivity, and more
# We'll manipulate the raster layer, National Land Cover Database (NLCD), to create separate raster layers important for amphibian presence.

# First, let's create matrices for the different layers we want.

# classifying the variables:
# 41, 42, 43  forest (frt): mixed, deciduous, evergreen
# 71          open land (old): shrub/scrub
# 90, 95      wetlands (wet): woody wetlands, emerging herbaceous
# 11          open water (owa): open water


mfrt <- matrix(c(0, 40, 0,
                 41, 43, 1, 
                 45, 99, 0), ncol = 3, byrow = T)
mold <- matrix(c(0, 70, 0,
                 71, 72, 1,
                 73, 99, 0), ncol = 3, byrow = T)
mwet <- matrix(c(0, 89, 0,
                 90, 95, 1,
                 96, 99, 0), ncol = 3, byrow = T)
mowa <- matrix(c(0, 10, 0,
                 11, 12, 1,
                 13, 99, 0), ncol = 3, byrow = T)

# classifying using terra::classify
frt <- classify(r["nlcd"], mfrt, right = NA)
old <- classify(r["nlcd"], mold, right = NA)
wet <- classify(r["nlcd"] mwet, right = NA)
owa <- classify(r["nlcd"], mowa, right = NA)

names(frt)<- "forest"
names(old)<- "open land"
names(wet)<- "wetlands"
names(owa)<- "open waters"

rlm <- c(frt, old, wet, owa)  # creating a raster stack 

plot(rlm)

# Let's crop the rlm objects to the study area
rlm.cp <- crop(rlm, study)
plot(rlm.cp)

# now, we're not calculating data for the entire state! Let's look at these values compared to our sites

df.buff <- st_buffer(df, dist = 500, joinStyle = "ROUND")
df.rlm <- terra::mask(x = rlm.cp, mask = df.buff, updatevalue = NA, touches = T)

tmap_mode("view")
tm_shape(df.rlm) +
  tm_raster(n=2,
            legend.show = T,
            labels = c(0,1),
            palette = c("white", "forestgreen"))+
  tm_shape(df)+
  tm_dots(size = 0.2)


# DEFINITIONS TO KNOW
# landscape metrics function at three levels: patch, class, and landscape. 
## patch = neighboring cells belonging to the same class, using the 8-neighbor rule 
## class = summary of all patches belonging to one class. Cna be a distribution of patch-level metrics of all patches (mean), or consider patches of the same class for calculations
## landscape = summarizes the whole landscape into one value

# CLASSES OF LANDSCAPE METRICS
# There are six different landscape metrics:
## Area and Edge Metrics: sizes of patches and classes and the amount of edge
## Shape Metrics: Shape of patch, mainly by area and perimeter
## Core Metrics: Area of patches that are not edge
## Aggregation Metrics: How patches of the same class are clumped or isolated from each other (i.e., the spatial configuration of the landscape)
## Diversity Metrics: Abundance and dominance/rareness of classes (landscape level only)
## Complexity metrics: How close to entropy the landscape pattern is (landscape level only)

# Let's take a look at our options:
View(list_lsm())

# First, pick a level. The level of landscape metric is dependent on your scientific question. 
## What is the goal of calculating landscape metrics for your study area? 
## What spatial scale are you trying to answer this question on?

# Also, note that many landscape metrics in the same class will (typically) be highly correlated to each other - pick landscape metrics from different classes 

# Finally, we need to think of buffer size around each point. How much does the landscape influence where a species/organism may (or may not be) present? 

# For this exercise, let's calculate some class-level metrics. 
## core area percentage of landscape (core area metric - lsm_c_cpland)
## patch density (aggregation metric - lsm_c_pd)
## aggregation (aggregation metric - lsm_c_ai) - for S&Gs


# creating a list of the variables we're interested in:
rlm.l <- list(rlm.cp["forest"], rlm.cp["wetland"], rlm.cp["open water"])

# running an lapply function for the landscape metrics we're interested in
psma.lm <- lapply(seq_along(rlm.l), function(x)
  sample_lsm(rlm.l[[x]], psma$geometry, level = "class", what = c("lsm_c_cpland", "lsm_c_pd", "lsm_c_ai"), shape = "square", size = 165, verbose = T, plot_id = psma$LabID))
  
# how is it calculating these landscape metrics? Via a *moving window*
# MOVING WINDOW: 
## A matrix that specifies the neighborhood and metric value for local neighborhood


# now, let's pivot the data so we can look at it

psma.lm <- lapply(seq_along(psma.lm), function(x)
  pivot_wider(psma.lm[[x]], names_from = metric, values_from = value))

View(psma.lm[[1]])

# note that there are some points that are duplicated - there are two class values per sample (1 and 0). What does this mean?

# If there are samples with two values, we care about the 1. Let's clean this up.
temp <- lapply(seq_along(psma.lm), function(x)
  psma.lm[[x]] %>% 
  group_by(plot_id) %>% 
  filter(n() == 2) %>% 
  filter(class == 1)
)



psma.lm <- lapply(seq_along(psma.lm), function(x)
  psma.lm[[x]] %>% 
  group_by(plot_id) %>% 
  filter(!n() == 2) %>% 
  bind_rows(temp[[x]]) %>% 
  arrange(plot_id) %>% 
  select(plot_id,
         cpland,
         pd,
         ai) 
)

# rename the variables so we know which variable is associated with what landscape
nms <- list("frt", "wet", "owa") # same order as the rlm vars

psma.lm <- lapply(seq_along(psma.lm), function(x)
  psma.lm[[x]] %>% 
    rename(setNames("cpland", paste("cpland", nms[[x]], sep = "_")),
           setNames("pd", paste("pd", nms[[x]], sep = "_")),
           setNames("ai", paste("ai", nms[[x]], sep = "_"))
           )
)

# check to make sure the names were transferred correctly
lapply(psma.lm, names)

# now, let's join this into one table

psma.dat <- psma.lm[[1]] %>% 
  right_join(psma.lm[[2]]) %>% 
  right_join(psma.lm[[3]])

View(psma.dat)

# now, for S&Gs, let's plot the relationship between two aggregate variables - the aggregation of patches (ai) and patch density (pd)

psma.dat %>% 
  ggplot(aes(x = ai_frt, y = pd_frt)) +
  geom_point()+
  geom_smooth(method = "lm") 

# let's look at the R-squared
summary(lm(ai_frt~pd_frt, data = psma.dat))

# note how correlated these two aggregation metrics are?
# Be careful when choosing landscape metrics - be sure to pick a metric from each category and not multiple from the same one! 
# These two do tell us slightly different things, but the base is the same


