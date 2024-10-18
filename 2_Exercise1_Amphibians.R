#_____________________________________________
# Exercise: Amphibians in SE Wyoming
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
                   "landscapemetrics",
                   "nlme",
                   "spmoran",
                   "graphics",
                   "Hmisc",
                   "corrplot",
                   "PerformanceAnalytics",
                   "spdep"))
# if rfUtilities is having problems:
remotes::install_github("jeffreyevans/rfUtilities")

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

# cool!
# So, our raster data is the entire state of Wyoming, whereas we are interested in shrinking our study area to capture just the points. Let's create a study area spatial file

# Let's see what the min and max is for our sites:
st_bbox(df)

# let's use the min/max of our points to create a polygon:
study <- st_as_sfc(st_bbox(df))

plot(study)
plot(df$geometry, add = T)

# Note that the bounding box is where the min/max points are. So, we'll increase the study area by adding a 1000m buffer around the bounding box
# This buffer ensures we have our data within the study area
# Creating a buffer also helps with preventing "edge effects" when we work with raster data - we don't want the boundary of our study area to be the reason why our analyses is biased

buff <- st_buffer(study, dist = 1000)

# cool! We have a study area! Let's crop our earlier rasters to our buffered study area:
rc <- terra::crop(r, buff)

# Now that our rasters are cropped to size, let's calculate a few more variables, specifically from our elevation layer (dem)
# slope (terra)
slope <- terrain(rc["dem"], v = "slope", neighbors = 8, unit = "degrees")

# aspect (terra)
asp <- terrain(rc["dem"], v = "aspect", neighbors = 8, unit = "degrees")

# we can also calculate HLI and SRR using the spatialEco package! However, for time's sake, I included those rasters from the beginning :)

# creating a stack of the rasters that have same extent, origin, crs, resolution
rs <- c(rc["dem"], rc["ffp"], rc["gsp"], rc["srr"], rc["hli"], slope, asp)

# Let's focus on one species for time's sake - boreal chorus frog (Pseudacris maculata, PSMA)
psma <- df %>% 
  filter(Species == "PSMA")

# let's also make sure we're getting one presence of PSMA per cell (this prevents over-inflating results and spatial autocorrelation)
# one way I've found success is to rasterize the points, then vectorizing again.
# creating a coordinates column 
coords <- st_coordinates(psma)

# rasterizing points to a reference raster, nlcd:
dat.r <- rasterize(cbind(coords[,1], coords[,2]), rc["nlcd"], value = 1, fun = mean)

# converting rasters into points
dat.rp <- as.points(dat.r)

# cleaning up the psma file for further processing
dat.psma <- st_as_sf(dat.rp) %>% 
  mutate(PointID = 1:nrow(dat.rp),
         Present = 1) %>% 
  select(-mean)

# we'll use the nlcd layer again from the `r` raster stack in step 2b.

#_____________________________________________
#### 2. Prepping Spatial Data for Analysis ####
#_____________________________________________
##### 2a. Creating Pseudoabsences #####
#_____________________________________________
# If I were to analyze the PSMA data as-is, the analyses assumes (a) perfect detection of PSMA and (b) I located all possible locations for PSMA in the study area. 
# I am good at finding frogs, but I'm not that good!
# So, we will prepare our data with that in mind!
# let's create pseudoabsences! 

rp <- st_sample(study, size = 2*nrow(dat.psma), type = "random") %>% 
  st_as_sf(data.frame(PointID = 1:(2*nrow(dat.psma)), Present = 0)) %>% 
  rename(geometry = x)

# let's add the random points to the presence points
dat.psma <- dat.psma %>% 
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
  mutate(Present = as.factor(Present),
         AllID = 1:nrow(dat.psma)) %>% 
  select(-ID) 

#_____________________________________________
##### 2c. Calculating Landscape Metrics #####
#_____________________________________________

# Landscape Metrics calculate the composition and configuration of different land uses within a region
# LMs are important for assessing habitat fragmentation/degradation, connectivity, and more
# We'll manipulate a raster layer, the National Land Cover Database (NLCD), to create separate raster layers important for amphibian presence.

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
frt <- classify(rc["nlcd"], mfrt, right = NA)
old <- classify(rc["nlcd"], mold, right = NA)
wet <- classify(rc["nlcd"], mwet, right = NA)
owa <- classify(rc["nlcd"], mowa, right = NA)

names(frt)<- "forest"
names(old)<- "open land"
names(wet)<- "wetlands"
names(owa)<- "open waters"

rlm <- c(frt, old, wet, owa)  # creating a raster stack 


# Let's look at these landuse values compared to our sites
df.buff <- st_buffer(df, dist = 500, joinStyle = "ROUND")
df.rlm <- terra::mask(x = rlm, mask = df.buff, updatevalue = NA, touches = T)

# we'll use tmap to observe these rasters
tmap_mode("view")
tm_shape(df.rlm) +
  tm_raster(n=2,
            legend.show = T,
            labels = c(0,1),
            palette = c("white", "forestgreen"))+
  tm_shape(df)+
  tm_dots(size = 0.2)

# Before we calculate the landscape metrics, we need to know the level and class of landscape metric
# LEVELS OF LANDSCAPE METRICS
## Think of levels like spatial scale - how much of the landscape do we want in our analysis? Does our species interact on a small scale or large scale?
## landscape metrics function at three levels: patch, class, and landscape. 
### patch = neighboring cells belonging to the same class, using the 8-neighbor rule 
### class = summary of all patches belonging to one class. Cna be a distribution of patch-level metrics of all patches (mean), or consider patches of the same class for calculations
### landscape = summarizes the whole landscape into one value

# CLASSES OF LANDSCAPE METRICS
## There are six different landscape metrics:
### Area and Edge Metrics: sizes of patches and classes and the amount of edge
### Shape Metrics: Shape of patch, mainly by area and perimeter
### Core Metrics: Area of patches that are not edge
### Aggregation Metrics: How patches of the same class are clumped or isolated from each other (i.e., the spatial configuration of the landscape)
## Diversity Metrics: Abundance and dominance/rareness of classes (landscape level only)
### Complexity metrics: How close to entropy the landscape pattern is (landscape level only)

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
## fractal dimension index (shape metric - lsm_c_frac_mn)
## aggregation (aggregation metric - lsm_c_ai) - for S&Gs


# creating a list of the variables we're interested in.
rlm.l <- list(rlm["forest"], rlm["open land"], rlm["wetland"], rlm["open water"])


# running an lapply function for the landscape metrics we're interested in
psma.lm <- lapply(seq_along(rlm.l), function(x)
  sample_lsm(rlm.l[[x]], dat.psma$geometry, level = "class", what = c("lsm_c_cpland", "lsm_c_pd", "lsm_c_frac_mn", "lsm_c_ai"), shape = "square", size = 165, verbose = T, plot_id = dat.psma$AllID))

# Let's see what this looks like
View(psma.lm[[1]])


# In this format, we can't do much. So, we'll pivot it to make more sense:
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
         frac_mn,
         ai) 
)


# rename the variables so we know which variable is associated with what landscape
nms <- list("frt","old", "wet", "owa") # same order as the rlm vars

psma.lm <- lapply(seq_along(psma.lm), function(x)
  psma.lm[[x]] %>% 
    rename(setNames("cpland", paste("cpland", nms[[x]], sep = "_")),
           setNames("pd", paste("pd", nms[[x]], sep = "_")),
           setNames("frac_mn", paste("fracc", nms[[x]], sep = "_")),
           setNames("ai", paste("ai", nms[[x]], sep = "_"))
    )
)

# check to make sure the names were transferred correctly
lapply(psma.lm, names)

# now, let's join the LMs into one table

psma.lm.all <- psma.lm[[1]] %>% 
  right_join(psma.lm[[2]]) %>% 
  right_join(psma.lm[[3]]) %>% 
  right_join(psma.lm[[4]])

psma.lm.all <- psma.lm.all %>% 
  rename(AllID = plot_id)

# now, for S&Gs, let's plot the relationship between two aggregate variables - the aggregation of patches (ai) and patch density (pd)

psma.lm.all %>% 
  ggplot(aes(x = ai_frt, y = pd_frt)) +
  geom_point()+
  geom_smooth(method = "lm") 

# let's look at the R-squared
summary(lm(ai_frt~pd_frt, data = psma.dat))

# note how correlated these two aggregation metrics are?
# Be careful when choosing landscape metrics - be sure to pick a metric from each category and not multiple from the same one! 
# These two do tell us slightly different things, but the base is the same

# Let's join these LMs to the main database
dat.psma <- dat.psma %>% 
  right_join(psma.lm.all)

# let's check for NAs
which(is.na(dat.psma), arr.ind = T)

# we have NAs in some rows/columns. Let's view these points and see what's up
na.dat <- dat.psma %>% 
  filter(if_any(everything(), is.na))

View(na.dat)

# what do these columns have in common?
# We can fortunately replace these values with a 0
dat.psma[is.na(dat.psma)] <- 0

# reordering variables - will be easier to call our variables this way
names(dat.psma)
dat.psma <- dat.psma %>% 
  relocate(AllID, .before = Present)

# Now, let's do some analyses!

#_____________________________________________
##### 3. Analyzing Spatial Data #####
#_____________________________________________
#_____________________________________________
##### 3a. Testing for Independence #####
#_____________________________________________
# We have the following:
## presence/pseudoabsence points
## covariates/predictor variables/independent variables***
## First, let's check to see how correlated our variables are by running a correlation matrix

c <- cor(st_drop_geometry(dat.psma[,-c(1:3)]), use = "complete.obs")

# that's gross. Let's filter out which variables are highly correlated with each other:

c[!lower.tri(c)] <- NA #removes diagonal and redundant values

cdf <- data.frame(c) %>% 
  rownames_to_column() %>% 
  gather(key = "variable", value="correlation", -rowname) %>% 
  filter(abs(correlation) > 0.7) # can change to be more or less conservative

View(cdf)

# now we can filter out which variables are redundant to each other!

# Below are a few ways to plot a correlation matrix that also show this relationshp
corrplot((cor(st_drop_geometry(dat.psma[,-c(1:3)]), method = "pearson")))

res <- rcorr(as.matrix(st_drop_geometry(dat.psma[,-c(1:3)])))

chart.Correlation(st_drop_geometry(dat.psma[,-c(1:3)]), histogram=TRUE, pch=19)
# last one is pretty ugly, but it works well when you're not dealing with many variables (maybe 5 - 6 vars)

# In general, if the absolute value of the correlation coefficient between variables =  0.7, then the two variables cannot be in the same model. 
# Pick one of the two variables, then move on.


# we can also use the `multi.collinearity` function in spatialEco to detect highly correlated variables:

# first, create a data frame of just the variables:
m <- dat.psma[,4:ncol(dat.psma)] 
m <- st_drop_geometry(m)
set.seed(425)

# inputing the data frame into the multi-collinearity function
mcl <- multi.collinear(m, perm = T, leave.out = T, n = 1000, p=0.05, na.rm = T)

View(mcl)

# testing for multi-collinearity is great because it assesses which variables are highly collinear with multiple variables at one time. The higher the frequency, the more times the variable was associated with another.

# To filter out our variables, we can remove them in clumps (for lack of a better word). We see a grouping of frequency = 700+, so let's filter those variables out first

m <- m[,-which(names(m) %in% mcl[mcl$frequency > 700,]$variables)]

# running again
mcl <- multi.collinear(m, perm = T, leave.out = T, n = 1000, p=0.05, na.rm = T)

View(mcl)

# Now, we have a grouping of frequency =00+
m <- m[,-which(names(m) %in% mcl[mcl$frequency > 400,]$variables)]

mcl <- multi.collinear(m, perm = T, leave.out = T, n = 1000, p=0.05, na.rm = T)

# We're at frequency = 0 for all variables! Now, we can continue!
vars <- mcl$variables

df <- dat.psma %>% 
  select(AllID, Present, all_of(vars))

#_____________________________________________
##### 3b. Linear Regression Models #####
#_____________________________________________
# Now that we have tested for variable independence, let's run a couple of models!
# We'll use a generalized linear model (GLM) because we're working with presence/absence data (binary). Linear models (LM) in R are not great at dealing with presence/absence data

# We have several variables, so we can hypothesize which variables influence PSMA presence on the landscape. 
names(psma.vars)

# water availability model
mod.wet <- glm(formula = Present ~ pd_wet, data = df %>% st_drop_geometry(), family = binomial(link="logit"))

# topography model
mod.top <- glm(formula = Present ~ dem + srr + slope + aspect, data = df %>% st_drop_geometry(), family = binomial(link="logit"))

# create your own model! What do you think?


# Now, let's get coefficients
# more parsimonious = less complex = better for models
# let's use AIC! Remember, lower is better

summary(mod.wet)
summary(mod.top)


# Let's do some model validation, specifically using Area under the Receiving Operator Curve (AUC and ROC)
# This model validation method tests the sensitivity and specificity (true positive vs false positive rate) of our model 
(wetAUC <- auc(dat.psma$Present, mod.wet$fitted.values))
(topAUC <- auc(dat.psma$Present, mod.top$fitted.values))


# our topographic model is the best one based on the AUC
# let's get the ROC plot to visualize the AUC
pred <- prediction(fitted(mod.top), df$Present)
perf <- performance(pred, measure="tpr", x.measure="fpr")
plot(perf, col=rainbow(10))
abline(coef=c(0,1))

# It's fine and dandy to create a model, but being able to predict back to the landscape is what we want. WHERE should we expect to have a good chance of PSMA occurring?

topr <- c(rs["dem"], rs["srr"], rs["slope"], rs["aspect"])

pred<- predict(topr, mod.top, type='response', progress="text")

ggplot() +
  layer_spatial(pred) +
  geom_sf(data = df, aes(group=as.factor(Present), color = as.factor(Present)))

# Note: if your model has landscape metrics, you need to create a raster of the landscape metric itself. You can use a moving window method to do such. It's very computationally intensive, so we're not going to do that for this exercise.

# Now we may want to choose some threshold: good for PSMA or not. 
# We'll create a probabilities data frame for where PSMA are found, along with a data frame of inputs vs fitted values
foundPSMA<-data.frame(obs = df$Present, fitted = mod.top$fitted.values) %>% 
  group_by(obs) %>% 
  summarise(Min = min(fitted), 
            Max = max(fitted))

# 0.3 (or 30% probability PSMA at spot) selected as threshold, so create a table we can use to reclassify the raster
rcls <- data.frame(from = c(0, .3), to = c(.3, 1), becomes = c(0, 1))

# Reclassify
predClass<- classify(pred, rcl = as.matrix(rcls))

ggplot() +
  layer_spatial(predClass) +
  geom_sf(data = df, aes(group=as.factor(Present), color = as.factor(Present)))

#_____________________________________________
##### 3c. Random Forest Models #####
#_____________________________________________
# Random Forest is a type of classification and regression tree (CART) model. 
# We can use this to determine the best set of variables that explain a species' presence on the landscape (versus our hypothesis method earlier).
# Let's try it out!

xdata <- df[,3:ncol(df)]
ydata <- df %>% 
  select(Present)

xdata <- st_drop_geometry(xdata)
ydata <- st_drop_geometry(ydata)

ydata <- as.factor(ydata$Present)

msel <- rf.modelSel(xdata = xdata, ydata = ydata, imp.scale="se", ntree = 501, r = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), final.model=T, seed = 425, replacement=F)

# look at this one first
msel$test



# this tells us all the models tested at the different importance percentiles (r in rf.modelSel)
# We are looking for (a) the lowest out of bag (OOB) error, lowest number of parameters, and lowest class error


# Let's continue with  a Random Forest model using "Selected variables:" and I prefer to use an odd number of trees to avoid ties.
modRf<- randomForest(x = df %>% 
                       st_drop_geometry() %>% 
                       dplyr::select(msel$selvars), y = df$Present, ntree = 501, importance = TRUE)

# view model summary
modRf

# What is the OOB error? Is that a high value or no? 
# What is the confusion matrix? What is that telling you?

# Let's continue and see the important factors
varImpPlot(modRf)

# Conduct internal prediction response
pred2<- predict(modRf, df %>% 
                  st_drop_geometry(), type ="response")


# internal prob response
obsProb<- as.data.frame(predict(modRf, df %>% 
                                  st_drop_geometry(), type = "prob"))
head(obsProb)

# make data frame with pred and prob
mObsPred<- data.frame(Observed = as.integer(as.character(df$Present)),
                      PRED = as.integer(as.character(pred2)),
                      Prob1 = obsProb[, 2],
                      Prob0 = obsProb[, 1])

# count correct preds
mop <-(mObsPred$Observed == mObsPred$PRED)
table(mop)

# validation rate
mpcc <- (length(mop[mop=="TRUE"])/length(mop))*100
mpcc

# internal auc. This isn't real, we need a training and a validation set. RF is too good when working with smaller datasets
intAUC<- auc(df$Present, obsProb[, 2])
intAUC

prdctnRf <- prediction(mObsPred$PRED, mObsPred$Observed)
perfRf <- performance(prdctnRf, measure = "tpr", x.measure = "fpr")
plot(perfRf, col = rainbow(10))
abline(coef = c(0, 1))


