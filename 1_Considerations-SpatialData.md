# Considerations for Spatial Data 
## Introduction
Before jumping off the deep end into the world of spatial data, having a foundation in both **Remote Sensing** and **Geographic Information Systems** (GIS) are important.

The information provided below **should not** be treated as a substitute for Remote Sensing/GIS courses. What I have provided in this GitHub/Workshop is introductory information for the lay viewer who might use spatial data in their work.

If your scientific question relies heavily on spatial data in any way, please take Remote Sensing/GIS courses. 

Also, while they are typically lumped together and there is information overlap, Remote Sensing and GIS are two different fields of study. Taking courses or certificates in both Remote Sensing and GIS will ensure you are receiving a holistic understanding of (a) spatial data sources, (b) understanding and interpreting spatial data, (c) manipulating and creating products from spatial data, and (d) analyzing spatial data.

## Behind the Scenes of Spatial Data: Terms to Understand
### *What is projection?*
* Spatial projection is the calculation used to flatten a 3D Earth onto a 2D plane
* There are different ways to mathematically calculate projections, and you are probably familiar with a few! Below are some projections.

<p>
<img src="https://github.com/user-attachments/assets/04a0db09-d01e-4c67-95d8-8344e557d8cb" width="25%" fig.show="hold">
<img src="https://github.com/user-attachments/assets/92cbde69-8d3e-45bf-8688-f91a527a403a" width="25%" fig.show="hold">
<img src="https://github.com/user-attachments/assets/6bebe2e4-67dc-450b-8c9c-3e4c786be2ed" width="25%" fig.show="hold">
<img src="https://github.com/user-attachments/assets/58e02b67-a9de-4374-a7a1-89bc3bcea900" width="20%" fig.show="hold">
  
  <h6><em> A selection of various projections. From Left to Right: Lambert Conformal Conic, Albers Equal Area, Wagner VI, and Web Map Mercator. Images sourced from Wikimedia Commons</em> </h6>
</p>

* Many projections exist because projecting Earth onto a 2D surface is not perfect - one projection may match one location perfectly, but not in another
* Knowing what projection you are using is important for ensuring your data (both raster and vector) align with each other

### *What are datums?*
* Projection is influenced by the **geodetic reference datum** (AKA datum), or the reference frame that positions locations on Earth via a set of coordinates. 
* Datums are calculated using a **reference ellipsoid**, or a simplified version of the Earth's shape whose equatorial axis is slightly longer than its polar axis.
* This ellipsoid is a good, but not perfect, estimate of the **geoid**, or the shape of Earth given its gravitational forces and rotation (not accounting for wind or tides). 
* Geoids are irregular in shape, whereas ellipsoids are not. So, different ellipsoids, thus different datums, are calculated globally and locally to improve GPS accuracy
* Some global datums are: **WGS84** and **ITRF**
* Some local datums are: **NAD83** and **NAD27** (North America), **ETRS89** (Europe)
* Also, datums have two components:
  * **Horizontal Datum**: the location of a spot. Uses coordinates.
  * **Vertical datum**: reference to vertical positions, like terrain, elevations, water level, etc.
* What projection, and ultimately what datum, you use are based on your research question and study area size.

### *What is Remote Sensing?*
* **Remote Sensing** is the process of capturing the landscape in a series of images that, typically, have data associated with them.
* Satellites, airplanes, and even drones can capture remotely sensed images
* Typically, imagery are captured in different **wavelengths**, or different bands in the light spectrum
* Different remotely sensed technologies have different sensors, thus can capture different scenes.
* **For example**, a drone might capture four bands (Red, Green, Blue, Near-Infrared), whereas NASA's Landsat 9 can capture 9 different bands, including thermal infrared
* This is important because we can create different imagery composites by manipulating the wavelengths we're seeing:

 <p align="center">
  <img src="https://github.com/user-attachments/assets/79e173ad-cc20-4c29-8908-57bdcf4e918d" width="75%">
<font align="center">
   <h6> 
  <em> Old Hickory Lake, TN, in two different spectrums. Left: Visible spectrum; Right: Near Infra-Red Spectrum (V. Yuksel, Wikimedia Common)</em> </font>
</h6>
</p>
 
* Manipulating different remotely sensed bands is how we can calculate different indexes, like Normalized Difference Vegetation Index (NDVI), calculate evapotranspiration on large scales, or even observe plant growth and stress.
* Thus, data sources are important to consider when answering your specific question!

### *What is Resolution?*
* **Resolution** is how much detail a remotely sensed image has
* There are four types of resolutions to consider. The broad definitions are:
  * **Spatial Resolution**: The size of the pixel within an image that represents an area
  * **Temporal Resolution**: The length of time a satellite takes to revisit the same observation area on its orbit. 
  * **Spectral Resolution**: The ability for a sensor to differentiate between wavelengths
  * **Radiometric Resolution**: The amount of information contained in each pixel.

You can visit [NASA's website here](https://www.earthdata.nasa.gov/learn/backgrounders/remote-sensing) to learn more about different resolutions

* While all resolutions are vital when deciding our raster data sources, if we are working in a spatial science project, I recommend thinking about *spatial* and *temporal* resolution first. 
  * For spatial resolution, the pixel size should represent (a) the study area, and (b) the species/dimension you are working with. **E.g.** a mule deer migrating thousands of miles vs a frog that migrates a maximum of 500m within an area will need different spatial resolutions to capture what influences movement for each species.
  * For temporal resolution, ensuring the satellite, plane, or drone collected the data during the right time of year, along with the frequency of that data collection, are important. Again, think about your scientific question: is the timing for when the data was collected important for your species/system? If so, ensure you are collecting spatial data during the correct time of year is important!
  

