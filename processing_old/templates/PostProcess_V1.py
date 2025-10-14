#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import pandas as pd
import time
import geopandas as gpd
import sys
import rasterio
from pathlib import Path

if (os.getcwd().split('/')[1] == 'Volumes') | (os.getcwd().split('/')[1] == 'Users'):
    path = '/Volumes/FortressL3/ANUGA/SWOT_Sim/'
    build_path = '/Users/Alchrist/Documents/Github/ANUGA/'
    code_path = '/Users/Alchrist/Documents/Github/ANUGA/processing/code/'
    templates_path = build_path + 'templates/'
    config_path = build_path + 'configure_files/'
    cnes_tools_path = build_path + 'swot-hydrology-toolbox/'#'/Volumes/FortressL3/ANUGA/GlobalDeltas/swot-hydrology-toolbox/'
    local_path = path + 'examples/'
else:
    path = '/projects/loac_hydro/alchrist/anuga/'
    #build_path = '/scratch_lg/loac_hydro/alchrist/anuga/build/'
    code_path = path + 'processing/code/'
    templates_path = path + 'templates/'
    config_path = path + 'configs/'
    cnes_tools_path = path + 'swot-hydrology-toolbox/'
    local_path = path + 'examples/'

sys.path.insert(1, code_path)
from delta_swot_simulations import (make_model_output_rasters,
                                    get_orbit_files ,
                                    run_swot_simulator,
                                    average_swot_pixel_clouds,
                                    convert_model_output_to_rasters)

from BYOM_extra import make_segments_for_swot
from BYOM_Utilities_V1 import (build_directory,
                               make_watermask,
                               make_channel_networks,
                               make_model_foundation)

AOI = '$AOI'
working_path,folders = build_directory(local_path, AOI)

modelpath = '$modelpath'
time_of_start = '$time_of_start'

parameters  = pd.read_csv(modelpath +'/' + 'config_%s_%s.csv' %(AOI,nowtime),dtype=str,delimiter=',')


EPSG = int(parameters['EPSG'][0])                                  # Coordinate system
elevation      = parameters['ElevationName'][0]                     # Name of the elevation file
elevationpath  = folders[4] + elevation
startdate      = (parameters['StartDate'][0])                    # Date to start simulation
enddate        = (parameters['EndDate'][0])                      # Date to start simulation
boundarytype1  = parameters['BoundaryType1'][0]                     # Boundary type for tidal boundary (See boundary types below - Bt sometimes crashes, Bout and Br are prefered)
boundarytype2  = parameters['BoundaryType2'][0]                     # Boundary type for remaining boundaires
discharge      = int(parameters['Discharge_cms'][0])                # River discharge in m3/s
meshname       = parameters['MeshName'][0]
scenario = '%s_%scms_%s_%s_%s%s%s%s' %(AOI,discharge,startdate,meshname,boundarytype2,boundarytype1,boundarytype1,boundarytype1)

EPSG = int(parameters['EPSG'][0])
mesh_scales     = parameters['Scale'].astype(int)                    # Maximum triangle area of meshes
base_scale      = min(mesh_scales) #50000                            # Maximum triangle area in ocean areas

startdate      = (parameters['StartDate'][0])                    # Date to start simulation
enddate        = (parameters['EndDate'][0])                      # Date to start simulation
xres           = (parameters['res'][0])                      # Date to start simulation

ulx  = parameters['ulx'][0]                                                 # ULX coordinate
lry  = parameters['lry'][0]                                                 # LRY coordinate
lrx  = parameters['lrx'][0]                                                 # LRX coordinate
uly  = parameters['uly'][0]

ref = rasterio.open('$dem')
ref.profile

superpixel_area = 10000
make_segments_for_swot(folders,AOI,ref,parameters,superpixel_area,False)

resampled_res = 30
time00 = time.time()

try:
    os.mkdir(folders[1])
except:''

AOI_file = gpd.read_file('%s%s_input.shp' %(folders[0],AOI))
AOI_file.crs = 'EPSG:4326'

cnespath = modelpath + '/SWOT_simulator/'
resampled_path = cnespath + str(resampled_res) +'/'
Path(cnespath).mkdir(parents=True, exist_ok=True)
print('CNES SWOT simulator files will be here: %s' %(cnespath))
Path(resampled_path).mkdir(parents=True, exist_ok=True)
Path('%s/orbits' %(cnespath)).mkdir(parents=True, exist_ok=True)

how_many_time_steps = 8
print('\n\n Convert SWW file to rasters for last 8 2-hr time steps \n\n')
convert_model_output_to_rasters(scenario,modelpath,nowtime,parameters,how_many_time_steps,resampled_res,resampled_path,False)

print('\n\n Making model output rasters for last 8 2-hr timesteps\n\n')
make_model_output_rasters(scenario,folders,nowtime,parameters,code_path,resampled_res,resampled_path,False)

print('\n\n Getting orbit files\n\n')
get_orbit_files(parameters,cnespath,AOI_file,cnes_tools_path,templates_path,os.path.isfile(cnespath+'/orbits/passplan.txt'))

print('\n\n Running SWOT simulator and stats analysis by segments\n\n')
run_swot_simulator(parameters,cnespath,folders,cnes_tools_path,templates_path,resampled_path,False)

print('\n\nAveraging simulated pixel clouds')
average_swot_pixel_clouds(parameters,cnespath,folders,modelpath,scenario,nowtime,superpixel_area,resampled_path,skip=False)
