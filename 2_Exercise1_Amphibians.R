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
# 



