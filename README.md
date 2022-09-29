# landspy
Landscape analysis with Python

## 1. Purpose

The landspy library provides some useful Python classes and functions for topographic analysis, drainage network extraction and computation of important geomorphic indices from Digital Elevation Models (DEMs). The main goal of landspy is to provide a simple way to compute geomorphic indices for those researchers who are not familiar with complex computational techniques. This library aims to pave the way for tectonic geomorphologists who do not have deep programming skills, but who can play a very relevant role in the interpretation of these analyses in terms of their relationship with geology, tectonics and associated hazards. 

This library is constantly growing, and we will include more analysis and functions in future versions. 

## 2. Installation

The library can be installed via pip via [pip][pip_link] or conda (we recommend to install via [conda-forge][conda_link]) on Linux, Mac, and Windows.

### pip installation
Install the package by typing the following command in a command terminal:

    pip install landspy

When installed via pip, it is recommended to have [GDAL][GDAL_pip_link] previously installed in our environment. It will be necessary to have libgdal and its development headers installed if pip is expected to do a source build because no wheel is available for your specified platform and Python version.

To install the latest development version via pip, see our [github project page][github_link].

### conda installation (recommended)
We can install landspy via [conda-forge][conda_forge_link].
Once you have Anaconda or Miniconda installed, you should be able to install GDAL with:

    conda install -c conda-forge landspy

### QGIS installation
When installing **landspy** to use inside QGIS, we will use the pip installation. As most of the requirements will be already available in the Python QGIS environment, follow these steps:
1. Open OSGeo4W shell (packed with QGIS in the start menu)
2. If running a QGIS version lower than 3.18 (not needed for later versions). Type:

    `py3_env`  
3. Install via pip by typing:

    `pip install landspy`

## 3. Requirements 
All dependencies should be installed along with landspy (See GDAL note for pip instalation).

- [GDAL](https://pypi.org/project/GDAL/)
- [NumPy >= 1.14.5](https://www.numpy.org)
- [SciPy >= 1.1.0](https://www.scipy.org/scipylib)
- [scikit-image >= 1.0.0](https://scikit-image.org/)
- [matplotlib >= 3.0.0](https://matplotlib.org/)

## 4. Citation

We are preparing a 
> Pérez-Peña et al.:
> Lansdpy, a open-source library for landscape analysis in Python and QGIS.
> [In preparation]

## 5. Tutorials and Examples

To get an overview of how **landspy** works, we offer some tutorials to perform some of the most common tasks that can be done with it. 

- [Extraction of a drainage network][tut1_link]
- [Calulation of Chi-maps and ksn values][tut2_link]
- [Plot single channel profiles][tut3_link]
- [Create channels from a polyline shapefile][tut4_link]
- [Create channels from a drainage basin][tut5_link]

The associated python scripts are provided in the `docs` folder.

### 6. Examples

#### Creation of a Chi Map from a Digital Elevation Model

This is an example of how to generate a Chi Map in vector format from a DEM.

```python
from landspy import DEM, Flow, Network
# Load the DEM and create the Flow and the Network
dem = DEM("data/jebja30.tif")
fd = Flow(dem)

# Create a Network object with a threshold of 1500 cells and reference m/n of 0.45
net = Network(fd, 1500, thetaref=0.45)

# Create a Chi Map in vector format for segments of 250 and 500 m
net.chiShapefile("data/chiMap_250.shp", 250)
net.chiShapefile("data/chiMap_500.shp", 500)
```

![ksn_values](https://user-images.githubusercontent.com/21242618/193000070-162ce11f-f729-49e4-b9cf-cbfe71461f62.jpg)

#### Analysis of different values of m/n for a basin

This is an example we analyze the best m/n value for Chi anlysis in a sample small basin

```python
from landspy import Grid, DEM, Flow, Network, BNetwork
import matplotlib.pyplot as plt 

# Load the DEM and create the Flow and the Network
dem = DEM("data/jebja30.tif")
fd = Flow(dem)

# Create a Network object with a threshold of 1500 cells and reference m/n of 0.45
net = Network(fd, 1500, thetaref=0.45)

# Load the Basins
basins = Grid("data/basins")

# Generate the BNetwork object for the basin with id=2
bnet = BNetwork(net, basins, bid =2 )

# Check different m/n values for chi analysis
mn_vals = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

fig = plt.figure(figsize=(17, 10))

for n, mn in enumerate(mn_vals):
    bnet.calculateChi(mn)
    ax = fig.add_subplot(2, 3, n+1)
    bnet.chiPlot(ax)
    ax.set_title("Chi plot (m/n = {})".format(mn))
    ax.set_xlabel("$\\chi$ (m)")
    ax.set_ylabel("Elevation (m)")
    
plt.tight_layout()
```

![mn_analysis](https://user-images.githubusercontent.com/21242618/193000280-b72ea0a5-8be8-47d6-a349-15a6053a8955.png)


## 7. Contact

You can contact me via <geolovic@gmail.com> <vperez@ugr.es>.

## 8. License

[MIT License][license_link] © 2022

[pip_link]: https://pypi.org/project/gstools
[conda_link]: https://anaconda.org/anaconda/repo
[GDAL_pip_link]: https://pypi.org/project/GDAL/
[conda_forge_link]: https://conda-forge.org/
[github_link]: https://github.com/geolovic/landspy
[tut1_link]: https://github.com/geolovic/landspy/tree/master/docs/landspy_tutorial_01.ipynb
[tut2_link]: https://github.com/geolovic/landspy/tree/master/docs/landspy_tutorial_02.ipynb
[tut3_link]: https://github.com/geolovic/landspy/tree/master/docs/landspy_tutorial_03.ipynb
[tut4_link]: https://github.com/geolovic/landspy/tree/master/docs/landspy_tutorial_04.ipynb
[tut5_link]: https://github.com/geolovic/landspy/tree/master/docs/landspy_tutorial_05.ipynb
[license_link]: https://github.com/geolovic/landspy/blob/master/LICENSE.txt