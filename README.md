# WorkShop-2024-R-SpatialScience
This repository is in collaboration with the R Ladies Workshop Series at the University of Wyoming. We'll review (a) what spatial science is, (b) considerations of spatial science, (c) best practices for collecting spatial data, and (d) learning how to model in R.

First, check out Spatial Science Considerations for a (very) brief overview of spatial data, including:
* projections and datums
* remote sensing
* resolutions

After reading the Considerations, check out Exercise 1

Also, please run the following code to ensure you have all packages:

```r{}
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
```

# References for this workshop
[Atkinson, P.M. & P.J. Curran.](https://openurl-ebsco-com.libproxy.uwyo.edu/EPDB%3Agcd%3A1%3A692449/detailv2?sid=ebsco%3Aplink%3Ascholar&id=ebsco%3Agcd%3A167357581&crl=c) 1997. Choosing an appropriate spatial resolution for remote sensing investigations. *Photogrammetric Engineering & Remote Sensing* 63(12):1345-1351

Dale, M.R.T and M.J. Fortin. 2014. Spatial Analysis: A Guide for Ecologists (Second Edition). Cambridge University Printing House.

[Lu, D. & Q. Weng.](https://www.tandfonline.com/doi/full/10.1080/01431160600746456) 2007. A survey of image classification methods and techniques for improving classification performance. *International Journal of Remote Sensing* 28(5):823-870

[Vermeer, G.J.O. 1999.](https://library.seg.org/doi/abs/10.1190/1.1444602) Factors affecting spatial resolution. *Geophysics* 64(3):942-953

With, K.A. 2019. Essentials of Landscape Ecology. Oxford University Press.











