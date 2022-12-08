#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANUGA Model
https://github.com/GeoscienceAustralia/anuga_core
Geoscience Australia, 2004-present

Code by Alexandra Christensen
NASA Jet Propulsion Laboratory
Decemeber 2020
Updated Version 3: 10-21-2021
Post-Process file
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
# Standard packages
import os
import pandas as pd
import time
import geopandas as gpd
import sys
from pathlib import Path

if os.getcwd().split('/')[1] == 'Volumes':
    path = '/Volumes/FortressL3/ANUGA/GlobalDeltas/'
    build_path = '/Users/Alchrist/Documents/Github/globalcoast/'
    code_path = '/Users/Alchrist/Documents/Github/globalcoast/code/'
    templates_path = build_path + 'templates/'
    config_path = build_path + 'configure_files/'
    cnes_tools_path = '/Volumes/FortressL3/ANUGA/GlobalDeltas/swot-hydrology-toolbox/'
    local_path = path + 'examples/'
else:
    path = '/scratch_lg/loac_hydro/alchrist/anuga/build/'
    build_path = '/scratch_lg/loac_hydro/alchrist/anuga/build/'
    code_path = '/scratch_lg/loac_hydro/alchrist/code/'
    templates_path = build_path + 'templates/'
    config_path = build_path + 'configure_files/'
    cnes_tools_path = '/scratch_lg/loac_hydro/alchrist/swot-hydrology-toolbox/'
    local_path = path + 'examples/'


sys.path.insert(1, code_path)
from delta_swot_simulations import make_model_output_rasters, get_orbit_files , run_swot_simulator,average_swot_pixel_clouds

modelpath = local_path + '$modelpath'
nowtime = '$nowtime'
parameters  = pd.read_csv(modelpath +'/' + '$config',dtype=str,delimiter=',')

delta = parameters['AOI'][0]                                       # Name of the AOI
deltapath = local_path + delta + '/'

folders = [deltapath + 'User_Defined_Files/',
           deltapath + 'tmp/',
           deltapath + 'Setup_Files/',
           deltapath + 'Meshes/',
           deltapath + 'Scenarios/',
           deltapath + 'Models/',
           deltapath + 'Setup_Files/Setup_SHP/',
           deltapath + 'Setup_Files/Setup_RST/',
           deltapath + 'Setup_Files/Setup_FIG/']
EPSG = int(parameters['EPSG'][0])                                  # Coordinate system
elevation      = parameters['ElevationName'][0]                     # Name of the elevation file
elevationpath  = folders[4] + elevation
startdate      = (parameters['StartDate'][0])                    # Date to start simulation
enddate        = (parameters['EndDate'][0])                      # Date to start simulation
boundarytype1  = parameters['BoundaryType1'][0]                     # Boundary type for tidal boundary (See boundary types below - Bt sometimes crashes, Bout and Br are prefered)
boundarytype2  = parameters['BoundaryType2'][0]                     # Boundary type for remaining boundaires
discharge      = int(parameters['Discharge_cms'][0])                # River discharge in m3/s
meshname       = parameters['MeshName'][0]
scenario = '%s_%scms_%s_%s_%s%s%s%s' %(delta,discharge,startdate,meshname,boundarytype2,boundarytype1,boundarytype1,boundarytype1)

time00 = time.time()

try:
    os.mkdir(folders[1])
except:''

AOI = gpd.read_file('%s%s_input.shp' %(folders[0],delta))
AOI.crs = 'EPSG:4326'


Path('%s/outputRST' %(local_path + modelpath)).mkdir(parents=True, exist_ok=True)
Path('%s/outputRST/abselev'%(local_path + modelpath)).mkdir(parents=True, exist_ok=True)
Path('%s/outputRST/depth'%(local_path + modelpath)).mkdir(parents=True, exist_ok=True)
Path('%s/outputRST/stage'%(local_path + modelpath)).mkdir(parents=True, exist_ok=True)
Path('%s/outputRST/waterelev'%(local_path + modelpath)).mkdir(parents=True, exist_ok=True)
Path('%s/outputRST/flooded'%(local_path + modelpath)).mkdir(parents=True, exist_ok=True)

cnespath = local_path + modelpath + '/SWOT_simulator'
Path(cnespath).mkdir(parents=True, exist_ok=True)
print('CNES SWOT simulator files will be here: %s' %(cnespath))

print('\n\n Making model output rasters for last 8 2-hr timesteps\n\n')
make_model_output_rasters(deltapath, code_path,scenario,local_path + modelpath,folders,cnespath,nowtime,parameters,templates_path,False)


print('\n\n Getting orbit files\n\n')
get_orbit_files(parameters,cnespath,AOI,cnes_tools_path, templates_path,os.path.isfile(cnespath+'/orbits/passplan.txt'))

print('\n\n Running SWOT simulator and stats analysis by segments\n\n')
run_swot_simulator(parameters,cnespath,folders,modelpath,scenario,cnes_tools_path,templates_path,nowtime,False)

# print('\n\n Averaging simulated pixel clouds')
# superpixels = gpd.read_file(folders[6] + '/gabon_oceanriver_bathy_superpixels_10000m2_1_0.95_felz_SWOT.shp')
# superpixel_area = 10000
# average_swot_pixel_clouds(parameters,cnespath,folders,modelpath,scenario,nowtime,superpixels,superpixel_area,skip=False)
