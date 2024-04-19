# SUDEM

Shielding Under a Digital Elevation Model (SUDEM) is a Matlab script that calculates the attenuation of cosmogenic radiation under an irregular surface.

The script calculates Shielding Factors for the different production mechanisms (spallation and muons) in samples located under the ground from a Digital Elvation Model and the location of the samples.

SUDEM is more accurate than considering a flat surface above the samples. The methodology is inspired by the method described in [Balco (2014)](https://doi.org/10.1016/j.quageo.2013.12.002). Although SUDEM is faster and easier to use, it is less precise when accounting for cosmic rays trajectories crossing the surface 3 or more times. Therefore, if samples are located under 'holes' (e.g. bedrock surface under a large erratic) Balco's method will produce more accurate results.

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
