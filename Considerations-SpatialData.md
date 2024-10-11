# Considerations for Spatial Data 
## Introduction
Before jumping off the deep end into the world of spatial data, having a foundation in both **Remote Sensing** and **Geographic Information Systems** (GIS) are important.

The information provided below **should not** be treated as a substitute for Remote Sensing/GIS courses. What I have provided in this GitHub/Workshop is introductory information for the lay viewer who might use spatial data in their work.

If your scientific question relies heavily on spatial data in any way, please take Remote Sensing/GIS courses. 

Also, while they are typically lumped together, Remote Sensing and GIS are two very different fields. Taking courses or certificates in both Remote Sensing and GIS will ensure you are receiving a holistic understanding of (a) spatial data sources, (b) understanding and interpreting spatial data, (c) manipulating and creating products from spatial data, and (d) analyzing spatial data

#### *What is projection?*
* Spatial projection is the calculation used to flatten a 3D Earth onto a 2D plane
* There are different ways to mathematically calculate projections, and you are probably familiar with a few! Below are some projections.

```{r out.width="25%", echo=F, message=F, warnings=F, fig.show="hold"}

knitr::include_graphics("C:/Users/melan/OneDrive/Desktop/drive/rworkshop_mlt/images/Lambert_conformal_conic_projection_SW.jpg")

knitr::include_graphics("C:/Users/melan/OneDrive/Desktop/drive/rworkshop_mlt/images/Albers_projection_SW.jpg")

knitr::include_graphics("C:/Users/melan/OneDrive/Desktop/drive/rworkshop_mlt/images/Wagner_VI_projection_SW.jpg")

knitr::include_graphics("C:/Users/melan/OneDrive/Desktop/drive/rworkshop_mlt/images/Web_maps_Mercator_projection_SW.jpg")

```
