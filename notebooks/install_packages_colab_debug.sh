#!/bin/bash


cd /content/drive/MyDrive/installations
echo "(1) Install pip packages to /content/drive/MyDrive/installations"
echo "nose mpi4py triangle dill Pmw pymetis mpi4py pyproj gdal geemap cmocean geopandas fiona pygeos rasterio rasterstats scikit-fmm rtree pyTMD Orinoco"
pip install -q --target=/content/drive/MyDrive/installations nose rtree shapely numpy mpi4py pyproj triangle dill Pmw pymetis mpi4py geemap cmocean geopandas fiona pygeos rasterio rasterstats scikit-fmm scipy backports.zoneinfo > /dev/null 2>&1

# echo "(2) Install pyTMD"
# pip install -q --target=/content/drive/MyDrive/installations git+https://github.com/tsutterley/pyTMD.git

wget -c https://repo.anaconda.com/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh
chmod +x Miniconda3-4.5.4-Linux-x86_64.sh
bash ./Miniconda3-4.5.4-Linux-x86_64.sh -b -f -p /usr/local
conda install -q -y --prefix /usr/local python=3.6 gdal


echo "(3) Install gdal"
apt-get -q -y install --target=/content/drive/MyDrive/installations python-gdal gdal-bin python3-gdal > /dev/null 2>&1
apt-get -y install --target=/content/drive/MyDrive/installations libgdal-dev
echo "(4) Install netcdf4"
apt-get -q -y install --target=/content/drive/MyDrive/installations python-netcdf4  > /dev/null 2>&1

echo "(5) Install OpenBlas"
apt-get -y install --target=/content/drive/MyDrive/installations libopenblas-dev > /dev/null 2>&1

echo "(6) Install Orinoco"
pip install --target=/content/drive/MyDrive/installations git+https://github.com/simard-landscape-lab/orinoco.git > /dev/null 2>&1

echo "(7) Download anuga_core github repository"
echo "htps://github.com/GeoscienceAustralia/anuga_core"
git clone https://github.com/GeoscienceAustralia/anuga_core.git  #> /dev/null 2>&1

echo "(8) Install anuga"
cd anuga_core
python setup.py --quiet build  > /dev/null 2>&1
python setup.py --quiet install  > /dev/null 2>&1

cd ../

echo "(7) Completed"