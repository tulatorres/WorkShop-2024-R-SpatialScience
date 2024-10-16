# Considerations for Spatial Data 
## Introduction
Before jumping off the deep end into the world of spatial data, having a foundation in both **Remote Sensing** and **Geographic Information Systems** (GIS) are important.

The information provided below **should not** be treated as a substitute for Remote Sensing/GIS courses. What I have provided in this GitHub/Workshop is introductory information for the lay viewer who might use spatial data in their work.

If your scientific question relies heavily on spatial data in any way, please take Remote Sensing/GIS courses. 

Also, while they are typically lumped together and there is information overlap, Remote Sensing and GIS are two different fields of study. Taking courses or certificates in both Remote Sensing and GIS will ensure you are receiving a holistic understanding of (a) spatial data sources, (b) understanding and interpreting spatial data, (c) manipulating and creating products from spatial data, and (d) analyzing spatial data.

## Behind the Scenes of Spatial Data: Terms to Understand
### *What is projection?*
* Spatial projection is the calculation used to transform a 3D object onto a 2D plane
* In most cases, this is referencing how to transform Earth (3D globe) onto a flat surface (2D, map)
* There are different ways to mathematically calculate projections, and you are probably familiar with a few!
* Figure 1 below shows some projections.

<p>
<img src="https://github.com/user-attachments/assets/04a0db09-d01e-4c67-95d8-8344e557d8cb" width="25%" fig.show="hold">
<img src="https://github.com/user-attachments/assets/92cbde69-8d3e-45bf-8688-f91a527a403a" width="25%" fig.show="hold">
<img src="https://github.com/user-attachments/assets/c55b64eb-4cfa-4d8f-83cb-e8c7cd993d0c" width="25%" fig.show="hold">
<img src="https://github.com/user-attachments/assets/58e02b67-a9de-4374-a7a1-89bc3bcea900" width="20%" fig.show="hold">
 

  <h6><em> Figure 1: A selection of various projections. From Left to Right: Lambert Conformal Conic, Albers Equal Area, Azimuthal Equidistant, and Web Map Mercator. (Image Source: Wikimedia Commons) </em> </h6>
</p>

* Many projections exist because projecting a 3D surface onto a 2D surface is not perfect - distortions of different features happen, especially with a complex 3D surface as Earth!
* For example, check out Figure 1 and compare Greenland vs Africa 
  *  In mercator projections that we are most familiar with (far right), Greenland looks approximately the same size as the continent, Africa.
  *  In reality, Africa (~11.7 mil mi<sup>2</sup> is significantly larger than Greenland (0.8 mil mi<sup>2</sup>)
* The most common projection methods you will encounter are: cylindrical, conic, and planar/azimuthal (Figure 2)

 <p align="center">
  <img src="https://github.com/user-attachments/assets/c67b9915-c99e-4413-bc44-918d64496756" width="50%">
<font align="center">
   <h6><em>Figure 2: Calculating different projections using cylindrical, conic, and azimuthal methods.(Image Source: Wikimedia Commons)</em></h6></font>
</p>


* In addition, many types of projections try to preserve certain properties of the 3D surface. These can be categorized as:
  * **Conformal**: perserving shape on a local spatial scale
  * **Equal-Area**: preserving area measurements
  * **Equidistant**: preserving distance between points
  * **Azimuthal**: preserving direction between points
* Knowing what projection is being used is important for ensuring your data (both raster and vector) align with each other




### *What are datums?*
* **geodetic reference datum** (AKA datum) is the mathematical model that provides reference points that positions locations on Earth via a set of coordinates. 
* Datums are calculated using a **reference ellipsoid**, or a simplified version of the Earth's shape whose equatorial axis is slightly longer than its polar axis.
* This ellipsoid is a good, but not perfect, estimate of the **geoid** (Figure 3), or the shape of Earth given its gravitational forces and rotation (not accounting for wind or tides).
  
 <p align="center">
  <img src="https://github.com/user-attachments/assets/4faaef3c-d53e-49a5-9e04-950230ab55a8" width="35%">
<font align="center">
   <h6> 
  <em> Figure 3: An example of a geoid (red line; uneven oval) vs the ellipsoid (black line, even oval) (Image Source: Wikimedia Commons) </em> </font>
</h6>
</p>

* Geoids are irregular in shape, whereas ellipsoids are not. So, different ellipsoids, thus different datums, are calculated globally and locally to improve GPS accuracy
* Some global datums are: **WGS84** and **ITRF**
* Some local datums are: **NAD83** and **NAD27** (North America), **ETRS89** (Europe)
* Also, datums have two components:
  * **Horizontal Datum**: the location of a spot. Uses coordinates.
  * **Vertical Datum**: reference to vertical positions, like terrain, elevations, water level, etc.
* What projection, and ultimately what datum, you use are based on your research question and study area size.

### *What are reference systems?*
* A **spatial reference system** (SRS), also known as a **coordinate reference system** (CRS), is what defines the coordinate system used to locate an object or feature
* Without an SRS/CRS, translating our data onto a map is impossible.
* SRS need three things:
  * Specific coordinate system (like lat/long, UTM, etc.)
  * A specified datum
  * A specified projection
* SRS/CRS are associated with both vector data (discrete data: points, lines, polygons) and raster data (continuous data: cells/pixels)
* Ultimately, the reference system ties the datum and the projection together - this is how we can interpret spatial data coordinates!
  * Example: an SRS might have datum WGS84 with projection UTM zone 13
  * Ensuring all data is aligned is important for (a) visualization and (b) appropriate analyses!

## Understanding GIS and Remote Sensing
### *What is GIS?*
* **Geographic Information Systems** (GIS) is the computing technology (both hardware and software) that captures, stores, manages, manipulates, visualizes, and analyzes geographical and spatial (i.e., geospatial) data.
* GIS is used every day to visualize and assess spatial patterns and processes
* Predicting storm trajectories, plotting species movements, assessing land use, predicting sea level rise, and more incorporate GIS methods
* GIS incorporates both vector and raster data for mapping and analyses
* 
### *What is Remote Sensing?*
* **Remote Sensing** is the process of capturing, visualizing, and analyzing the landscape via imagery captured aerially
* Satellites, airplanes, and even drones can capture remotely sensed images
* These images are rasters - pixels/cells with data associated with it
* Typically, imagery are captured in different **wavelengths**, or different bands in the light spectrum
* Different remotely sensed technologies have different sensors, thus can capture different scenes.
* **For example**, a drone might capture four bands (Red, Green, Blue, Near-Infrared), whereas NASA's Landsat 9 can capture 9 different bands, including thermal infrared
* This is important because we can create different imagery composites by manipulating the wavelengths we're seeing:

 <p align="center">
  <img src="https://github.com/user-attachments/assets/79e173ad-cc20-4c29-8908-57bdcf4e918d" width="75%">
<font align="center">
   <h6> 
  <em> Figure 3: Old Hickory Lake, TN, in two different spectrums. Left: Visible spectrum; Right: Near Infra-Red Spectrum (Image Source: V. Yuksel, Wikimedia Commons)</em> </font>
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
 

  

