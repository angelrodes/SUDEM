# SUDEM

Shielding Under a Digital Elevation Model (SUDEM) is a Matlab script that calculates the attenuation of cosmogenic radiation under an irregular surface (self-shielding).

The script calculates Shielding Factors for the different production mechanisms (spallation and muons) in samples located under the ground from a Digital Elvation Model and the location of the samples.

The mehod is based in the interpolation of the given DEM to generate radial DEMs with log-distributed distances from the sample. Cosmic rays through the sample and each node of the radial DEM are simulated.

Shielding factors calculated using SUDEM are more accurate than considering a flat surface above the samples. The methodology is inspired by the code described in [Balco (2014)](https://doi.org/10.1016/j.quageo.2013.12.002). Although SUDEM is faster and easier to use, it is less precise when accounting for cosmic rays trajectories crossing the surface 3 or more times. Therefore, if samples are located under 'holes' (e.g. bedrock surface under a large erratic) Balco's method produces more accurate results.

[Angel Rodes](http://www.angelrodes.com), 2024.

## How to use SUDEM

SUDEM requires two files:

### Digital Elevation Model

The digital elevation model of the area should be exported or downloaded (e.g. 5-m Digital Models of the Terrain from [Digimap](https://digimap.edina.ac.uk/).[Digimap](https://digimap.edina.ac.uk/)) as a text ```.xyz``` file, renamed to ```dem.xyz``` and placed in the same folder as SUDEM script.

The ```dem.xyz``` file is a simple text-based format for storing elevation data in a three-dimensional coordinate system. Each line in the file represents a single point in space, with the X, Y, and Z coordinates separated by commas. Here's what each component typically represents:

-   X: The horizontal position in the coordinate system, usually representing longitude or easting.
-   Y: The vertical position in the coordinate system, usually representing latitude or northing.
-   Z: The elevation or height above a reference surface, usually sea level.

**X, Y, and Z must be expressed in metres** for SUDEM to calculate the attenuation of the cosmic radiation properly.

### Sample locations

A second file containing the sample locations is also required in the SUDEM folder. The ```samples.xyz``` file has the same format and units as ```dem.xyz```. However, the elevations (Z) in the ```samples.xyz``` file can be expressed as absolute (relative to sea level, as in ```dem.xyz```) or relative to the ```dem.xyz``` surface (using values of 0 or negative for Z). SUDEM will atomatically detect how the elevations are expressed as following:

- If any of the Z values is positive, all the samples are considered to be expressed as absolute elevations above sea level.
- If all the Z values are 0 or negative:
    - They will be considered as elevations relative to the ```dem.xyz``` surface if the minimum elevation of the ```dem.xyz``` surface is at least 100 mabove sea level.
    - Otherwise, a dialog will ask the user to choose if the Z values are realtive or absolute.

**Note that X, Y, and Z values must also be expressed in metres.** For depth profiles, remember to convert to negative metres!

## Other parameters

At the top of the Matlab script, the following parameters can be changed:

- Number of points for radial DEMs: Resolution of the interpolated digital elevation models. 100 azimuths and 300 distances are considered by default. More resolution, more precission, longer runs.
- ```min_distance```: Minimun horizontal distance for near vertical cosmic rays. 10 cm is considered by default. For very deep samples, consider reducing this number.
- ```randomized_models```: % Number of randomized model used to estimate uncertainties.
- Attenuationg lengths for cosgenic productions: 160, 850, 5000, and 500 g/cm2, according to the aproximations described in [Rodés, Á. (2021)](https://doi.org/10.3390/geosciences11090362).
- Density of the rock.
- ```display_3D_dem```: Choose if showing (```1```) or not (```0```) the DEM in a 3-D plot. ```0``` will save computing time.

## Otputs

The scipt will produce a figure for each sample. Each figure show the different surface transects radial from the sample location. As the distances are dirtibuted logarithmically, two plots are produced: one for local morphology (angles<20º), and one for the entire radial DEM.

![image](https://github.com/angelrodes/SUDEM/assets/53089531/13b52b60-11e4-4b87-9718-69ea99a6f7f8)
*Colour hue of the cross-section correspond to azimuths.*

A table showing all the shielding factors and uncertainties will be shown in the Command Window.

![image](https://github.com/angelrodes/SUDEM/assets/53089531/7f35e3d2-f131-4eda-8aa6-06d6c2e73c54)

## Under the hood

![SUDEM_story](https://github.com/angelrodes/SUDEM/assets/53089531/ba4bde1c-4cd8-4a9d-8a1a-99a992f81004)

### Radial DEM

SUDEM genrates radial grids centered in the sample locations. The nodes of the grid are ditributed along ```n_azimuths``` directions and ```n_distances``` distances spaced logartithmically, form a horizontal ```min_distance``` from the sample to the longest distance from the sample with an elevation higher than the sample in the original DEM. Node elevations are 2-D interpolated from the original DEM. These nodes are then used to simulate cosmic rays arriving to the sample.

### Shielding calcualtions

For each ray trajectory, the model calculates the distance (d) between the sample and the surface and the zenith angle (φ). The attenuation of the cosmic ray for this trajectory is calculated as e<sup>-d*ρ/Λ</sup>, where ρ is the density of the rock and Λ is the attenuation length of the cosmic radiation.

All the attenuations calculated are weighted by N<sub>ΔΦ</sub> · _sin_<sup>2.3</sup>(φ), where N<sub>ΔΦ</sub> is the number of cosmic ray for each of the 100 zenith angle groups (one every 0.9°), and _sin_<sup>2.3</sup>(φ) is the relative contribution of the cosmic radiation according to [Gosse & Phillips (2001)](https://doi.org/10.1016/S0277-3791(00)00171-2). The final shielding factor is calculated as 1 minus the weighted average attenuation.

In contrast to [Balco (2014)](https://doi.org/10.1016/j.quageo.2013.12.002)'s method, holes that make cosmic rays to cross the surface 3 or more times are not simulated, and their contributions are considered "in average" according to the different distances between the surface and the sample. For a smooth surface, we don't expect many close-to-vertical rays, the ones that contribute more, to cross the surface 3 or more times. Therefore, this approximation provides accurate shielding factors for most applications.

### Uncertainties

In order to estimate the uncertainty of this method, the radial grids are re-calulated with nodes corresponding to random interpolations between the initial nodes of the radial DEM. Thus, azimuths and distances of randomized grids will not be homogeneously distributied. The first random distance is randomized between 0 and ```min_distance```, so randomized grids have the same number of nodes as the original grid. The standard deviation of the the produced shielding factors is shown as uncertinty.
 
