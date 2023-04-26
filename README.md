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

Several datasets are utilized in this project. Those that cannot be accessed through a Python API must be downloaded manually by the user. Users must populate the BAM/ancillary folder with the files below for the scripts to work. They are too large to be provided through GitHub.

### GEBCO General Bathymetry Charts of the Ocean
> ***Download the files from here: https://www.bodc.ac.uk/data/open_download/gebco/gebco_2022_sub_ice_topo/geotiff/, unzip the file, and move to BAM/ancillary/gebco. There is a vrt file there that will only work with this version of GEBCO.
 - The 2022 global gridded dataset are subset within this project. The global database is large (8GB) and can be slow to download
 - Information https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/
 - Citation: GEBCO Compilation Group (2022) GEBCO_2022 Grid (doi:10.5285/e0f0bb80-ab44-2739-e053-6c86abc0289c)
 
### Global Mangrove Watch (BAM/ancillary/GMW)
> ***Download these v3_2016 file https://zenodo.org/record/6894273/files/gmw_v3_2016_vec.zip?download=1, unzip it, and move it to BAM/ancillary/GMW***
 - Mangrove habitat extent - we use 2016 shapefiles in this project. Version 3 is now available
 - All downloads are here: https://zenodo.org/record/6894273. 
 - Interactive map available here: https://www.globalmangrovewatch.org/
 - Citation: Bunting P, Rosenqvist A, Hilarides L, Lucas RM, Thomas N, Tadono T, Worthington TA, Spalding M, Murray NJ, Rebelo L-M. Global Mangrove Extent Change 1996–2020: Global Mangrove Watch Version 3.0. Remote Sensing. 2022; 14(15):3657. https://doi.org/10.3390/rs14153657 

### World Water Bodies
- Available from UCLA
- Download here: https://apps.gis.ucla.edu/geodata/hr/dataset/world_water_bodies/resource/a6b40af0-84cb-40ce-b1c5-b024527a6943

### Geoid Corrections
> ***Download this EGM2008 file https://vdatum.noaa.gov/download/data/vdatum_EGM2008.zip, unzip it, find the file "egm2008.gtx" and move it to BAM/ancillary/geoids***
- These gtx files are used to convert between different geoids
- Available from NOAA's VDatum (resource for converting between vertical datums)
- All gtx files are here: https://vdatum.noaa.gov/download.php#download_geoid

 <br></br>



## Tutorials: 
Tutorials can be run on your local computer or through Google Colab using the "Open in Colab" buttons below. The tutorials use the Komo Estuary in Gabon as an example. More tutorials will be uploaded soon. <br>

### Tutorial 0: Google Colab Introduction
Summary: This notebook will walk through the process of installing packages, connecting to your Google Drive, signing up for Google Earth Engine, and using Google Colab (a free, online platform for running Python Jupyter Notebooks)<br></br>
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/achri19/BAM/blob/main/notebooks/GoogleColab_Introduction.ipynb)<br></br>

### Tutorial 1: Build Model Files
Summary: This notebook will walk through steps to build a digital elevation for the study area using open-source/publicly available datasets.<br></br>
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/achri19/BAM/blob/main/notebooks/Build_Model_Files.ipynb)
<br></br>




## Citation:
Author: Alexandra Christensen
Affiliation: Jet Propulsion Laboratory, California Institute of Technology
Acknowledgement: The research was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (80NM0018D0004)

Copyright 2022, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.

This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons
