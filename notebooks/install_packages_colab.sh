#!/bin/bash


cd /content/drive/MyDrive/installations

echo "(1) Install pip packages to /content/drive/MyDrive/installations"
echo "nose mpi4py triangle dill Pmw pymetis mpi4py pyproj gdal geemap cmocean geopandas fiona pygeos rasterio rasterstats scikit-fmm rtree pyTMD Orinoco"
pip install --upgrade --target=/content/drive/MyDrive/installations nose rtree shapely numpy pandas  mpi4py triangle dill Pmw pymetis geemap cmocean geopandas netCDF4 pyproj fiona pygeos rasterio rasterstats scikit-fmm scipy  #> /dev/null 2>&1


echo "(2) Install gdal"
apt-get -y install gdal-bin python3-gdal #> /dev/null 2>&1

echo "(3) Install Orinoco"
pip install --target=/content/drive/MyDrive/installations git+https://github.com/simard-landscape-lab/orinoco.git > /dev/null 2>&1

echo "(4) Download anuga_core github repository"
echo "htps://github.com/GeoscienceAustralia/anuga_core"
git clone https://github.com/GeoscienceAustralia/anuga_core.git  #> /dev/null 2>&1

echo "(5) Install anuga"
cd anuga_core
python setup.py --quiet build  > /dev/null 2>&1
python setup.py --quiet install  > /dev/null 2>&1

cd ../

echo "(6) Completed"