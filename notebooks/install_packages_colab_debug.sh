
#!/bin/bash

mkdir -p /content/drive/MyDrive/installations
cd /content/drive/MyDrive/installations
echo "(1) Install pip packages to /content/drive/MyDrive/installations"
echo "nose mpi4py triangle dill Pmw pymetis mpi4py gdal netcdf4 pyproj geemap cmocean geopandas fiona pygeos rasterio rasterstats scikit-fmm rtree pyTMD Orinoco"
pip install --target=/content/drive/MyDrive/installations nose mpi4py pyproj triangle dill Pmw pymetis mpi4py netcdf4 geemap cmocean geopandas fiona pygeos rasterio rasterstats scikit-fmm rtree backports.zoneinfo pyTMD git+https://github.com/simard-landscape-lab/orinoco.git #> /dev/null 2>&1

echo "(2) Install gdal"
apt-get -y install --target=/content/drive/MyDrive/installations python-gdal gdal-bin  > /dev/null 2>&1

echo "(3) Install netcdf4"
apt-get -y install --target=/content/drive/MyDrive/installations python-netcdf4  > /dev/null 2>&1

echo "(3) Download anuga_core github repository"
echo "https://github.com/GeoscienceAustralia/anuga_core"
git clone https://github.com/GeoscienceAustralia/anuga_core.git  #> /dev/null 2>&1

echo "(4) Install anuga"

cd anuga_core
python setup.py build #--quiet build  > /dev/null 2>&1
python setup.py install #--quiet install  > /dev/null 2>&1

cd ../

echo "(7) Completed"
