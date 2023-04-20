#!/bin/sh
eval "$start"
conda activate swot_env2
conda info|grep "active environment"

export SWOT_HYDROLOGY_TOOLBOX=/users/alchrist/documents/tools/swot-hydrology-toolbox
export RIVEROBS=/users/alchrist/documents/tools/RiverObs

#cd '$command'#/Users/Alchrist/Documents/Tools/swot-hydrology-toolbox/test/guayas/guayas_1047cms_20210101_Meshes_river_10000m2_ocean_50000m2_land_20000m2_lake_10000m2_BoutBrBrBr_20211221/data2/
files='$files'
for f in files
do
  ogr2ogr -F netCDF '$command'
done


## pip uninstall shapely
## pip install shapely --no-binary :all:
