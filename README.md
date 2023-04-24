# BAM: Building ANUGA Models 

This repository contains Python scripts for building and running open-source hydrodynamic ANUGA models using publicly available data. More information about ANUGA is available here https://github.com/GeoscienceAustralia/anuga_core

## Installation

Use the BAM.yml file to set up a conda environment with requried packages:
```diff
- conda env create -f BAM.yml 
```

If the installation stalls at "Solving Environment", we recommend using Mamba (https://github.com/mamba-org/mamba) to build the environment. From base environment:

``` diff
- conda install mamba -n base -c conda-forge
- mamba env create -f BAM.yml
```

To activate the new environment:
``` diff
- conda activate BAM
```


## Ancillary Datasets

Several datasets are utilized in this project. Those that cannot be accessed through a Python API must be downloaded manually by the user.

GEBCO General Bathymetry Charts of the Ocean
 - The 2020 or 2022 global gridded dataset can be downloaded here (8 GB)
 - Information https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/
 - Download here: https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_sub_ice_topo/geotiff/
 - Citation: GEBCO Compilation Group (2022) GEBCO_2022 Grid (doi:10.5285/e0f0bb80-ab44-2739-e053-6c86abc0289c)
 
 

Author: Alexandra Christensen
Affiliation: Jet Propulsion Laboratory, California Institute of Technology
Acknowledgement: The research was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (80NM0018D0004)
 <br></br>



## Tutorials: 
Tutorial 0: Google Colab Introduction <br></br>
Summary: This notebook will walk through the process of installing packages, connecting to your Google Drive, signing up for Google Earth Engine, and using Google Colab (a free, online platform for running Python Jupyter Notebooks)<br></br>
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/achri19/BAM/blob/main/notebooks/GoogleColab_Introduction.ipynb)<br></br>

Tutorial 1: Build Model Files <br></br>
Summary: This notebook will walk through steps to build a digital elevation for the study area using open-source/publicly available datasets.<br></br>
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/achri19/BAM/blob/main/notebooks/Build_Model_Files.ipynb)
<br></br>



<br></br><br></br><br></br>
Copyright 2022, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.

This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons
