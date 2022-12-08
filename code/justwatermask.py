#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 07:28:56 2021

@author: Alchrist
This script is specifically for developing new SWORD channel make_SWORD_networks
It used a 10m water mask build in GEE from optical Sentinel-2 and/or radar Sentinel-1 images
It can be run locally or on Gattaca, but paths should match those below

"""

import sys
import os
import pandas as pd
import shutil
from datetime import datetime
from string import Template
import fnmatch
import geopandas as gpd
import platform
import psutil
import rasterio
from osgeo import gdal
from pathlib import Path
import numpy as np

if os.getcwd().split('/')[1] == 'Users':
    path = '/Users/alchrist/Documents/GitHub/ANUGA/'
    path_code = path + 'processing/code/'
    path_templates = path + 'processing/templates/'
    path_configs = path  + 'processing/configs/'
    cnes_tools_path = path + 'swot-hydroloy-toolbox/'
    path_ancillary = path + 'ancillary/'
    path_examples = '/Volumes/FortressL3/ANUGA/GlobalDeltas/examples/' #''/Volumes/FortressL3/ANUGA/GlobalDeltas/examples/'
else:
    path = '/projects/loac_hydro/alchrist/anuga/'
    path_code = path + 'code/'
    path_templates = path + 'templates/'
    path_configs = path + 'configs/'
    path_ancillary = path + 'ancillary/'
    cnes_tools_path = '/scratch_sm/loac_hydro/swot-hydrology-toolbox/'
    path_examples = path + '/examples/'

sys.path.insert(1,path_code)
""" These functions were developed for the BYOM tutorials, but are used for developing """

from BYOM_Utilities_V1 import (build_directory,
                               get_extent_parameters,
                               make_polygons,
                               make_channel_networks,
                               make_watermask,
                               more_opening)

print(f"Computer network name: {platform.node()}") #Computer network name
print(f"Machine type: {platform.machine()}") #Machine type
print(f"Processor type: {platform.processor()}") #Processor type
print(f"Platform type: {platform.platform()}") #Platform type
print(f"Operating system: {platform.system()}") #Operating system
print(f"Operating system release: {platform.release()}") #Operating system release
print(f"Operating system version: {platform.version()}") #Operating system version
print(f"Number of physical cores: {psutil.cpu_count(logical=False)}") #Physical cores
print(f"Number of logical cores: {psutil.cpu_count(logical=True)}") #Logical cores
print(f"Current CPU utilization: {psutil.cpu_percent(interval=1)}") #System-wide CPU utilization
print(f"Current per-CPU utilization: {psutil.cpu_percent(interval=1, percpu=True)}") #System-wide per-CPU utilization
print(f"Total RAM installed: {round(psutil.virtual_memory().total/1000000000, 2)} GB") #Total RAM
print(f"Available RAM: {round(psutil.virtual_memory().available/1000000000, 2)} GB") #Available RAM
print(f"Used RAM: {round(psutil.virtual_memory().used/1000000000, 2)} GB") #Used RAM
print(f"RAM usage: {psutil.virtual_memory().percent}%") #RAM usage

##############################################################################
######################### Get Configuration Parameters #######################
##############################################################################
# Resolution all output files - ANUGA models are always in UTM
res = 10

configfile = '//Users/Alchrist/Documents/Github/globalcoast/configs/Deltas_SetupParameters_MASTER_V2.csv' # #'Deltas_SetupParameters_MASTER_V2.csv'path_configs + 'need_opening.csv'#'
parameter_file = pd.read_csv(configfile)
AOIs = parameter_file['AOI']

i=0
#AOIs = ['guayas']

print('\n\n\n##############################################################################################')
print('#################################[Step 3A][Make_Watermask]####################################')
print('##############################################################################################\n')

##############################################################################
##################### Get best watermask from GEE ############################
##############################################################################
print('\n[Step 2][Setup_AOI_Files][Compiling Water Masks] .......\n')
#gee_ndwi_path = path + 'GEE_NDWI_watermask' #'/Volumes/GoogleDrive/My Drive/GEE_drive/'
# gee = [os.path.join(dirpath,f)
#     for dirpath,dirnames, files in os.walk(gee_ndwi_path)
#     for f in fnmatch.filter(files,'%s*' %(delta))]
#if len(gee)<=0:
xres = 10

hydropolys = gpd.read_file("%sucla/hydropolys_fix.shp" %(path_ancillary))
plds = gpd.read_file('%sSWOT_PLD.gdb/SWOT_PLD.shp' %(path_ancillary))
hydrolakes = gpd.read_file('/Volumes/FortressL3/ExternalData/HydroSheds/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_fixed.shp')
pld_crs = plds.crs
hydrolakes_crs = hydrolakes.crs
hydropolys_crs = hydropolys.crs
from shapely.geometry import Polygon, LineString

for AOI in AOIs[:]:
    skip = False
    print('################ ################# ##################')
    print('################ ################# ##################')
    print('################ ################# ##################')
    print('################ ANUGA Model Setup ##################')
    print('################ ################# ##################')
    print('################ ################# ##################')
    print('################ ################# ##################')


    AOI = AOI.lower()
    parameters = parameter_file[parameter_file['AOI']==AOI.capitalize()].reset_index(drop=True)

    ## Set up working directory and build sub folders if needed
    print('\n')
    print('Study area is ' + AOI)
    print('Resolution of this setup is %sm' %(res))
    working_path,folders = build_directory(path_examples, AOI)

    final_watermasks = [os.path.join(dirpath,f)
        for dirpath,dirnames, files in os.walk(folders[0])
        for f in fnmatch.filter(files,'*_finalwatermask*.tif')]

    print('Before building the model, we need 3 pieces: \n1) model domain (shapefile) \n2) approximate tide boundary (shapefile) \n3) watermask')
    #if len(final_watermasks) >0:
    print('You already have a chosen water mask:')
    print(final_watermasks)

    ref_10m,parameters = get_extent_parameters(path_ancillary,AOI,folders,res,parameters)

    print(parameters.iloc[0] )

    EPSG = parameters['EPSG'][0]
    ulx = parameters['ulx'][0]
    uly = parameters['uly'][0]
    lrx = parameters['lrx'][0]
    lry = parameters['lry'][0]

    extentpoly = gpd.read_file("%s%s_extent_%s.shp" %(folders[7],AOI,EPSG))

    model_domain = gpd.read_file('%s%s_modeldomain.shp' %(folders[7],AOI))
    AOI_extent = gpd.read_file('%s%s_extent_%s.shp' %(folders[7],AOI,EPSG))

    clean_with_landcover = False
    print('\n[Step 3A][Make_Watermask][Load Hydropolys] .......\n')
    extentpoly2 = extentpoly.buffer(xres*10)                                   # Create slightly larged extent for clipping in 4326 crs
    extentpoly2 =  extentpoly2.to_crs(hydropolys_crs)                       # Reproject to 4326 CRS
    hydropoly = gpd.clip(hydropolys,extentpoly2)                            # Clip hydropolys to model domain extent
    hydropoly = hydropoly.to_crs("EPSG:%s" %(EPSG))                         # Reproject to UTM
    hydropoly = gpd.overlay(hydropoly,extentpoly,how='intersection')        # Clip to extent in correct EPSG, ensure clipped area is correct
    hydropoly.to_file("%s%s_hydropolys.shp" %(folders[0],AOI))            # Save

    print('##################### Extracting SWOT PLD Lakes for the AOI')
    extentpoly2 = extentpoly.buffer(xres)                                   # Create slightly larged extent for clipping in 4326 crs
    extentpoly2 =  extentpoly2.to_crs(pld_crs)                       # Reproject to 4326 CRS
    pld = gpd.clip(plds,extentpoly2)                            # Clip HydroLakes to model domain extent
    pld = pld.to_crs("EPSG:%s" %(EPSG))                         # Save
    if len(pld) == 0:
        print('There are no lakes in the model domain')
        d = {'geometry': [Polygon([(0, 0), (0,0),(0,0)])]}
        pld = gpd.GeoDataFrame(d,crs="EPSG:%s" %(EPSG))
    else:
        pld['TYPE_2'] = 'Lake'                                            # Add label to identify lakes
        pld.geometry = pld.buffer(0)                               # Extra buffer to account low precision of hydrolakes
        pld = gpd.overlay(pld,extentpoly,how='intersection')        # Clip to extent in correct EPSG, ensure clipped area is correct
    pld.to_file("%s%s_SWOTPLD.shp" %(folders[0],AOI))

    print('##################### Extracting HydroLAKES for the AOI')
    hydrolake = gpd.clip(hydrolakes,extentpoly2)                            # Clip HydroLakes to model domain extent
    hydrolake = hydrolake.to_crs("EPSG:%s" %(EPSG))                         # Save
    if len(hydrolake) == 0:
        print('There are no lakes in the model domain')
        d = {'geometry': [Polygon([(0, 0), (0,0),(0,0)])]}
        hydrolake = gpd.GeoDataFrame(d,crs="EPSG:%s" %(EPSG))
    else:
        hydrolake['TYPE_2'] = 'Lake'                                            # Add label to identify lakes
        hydrolake.geometry = hydrolake.buffer(50)                               # Extra buffer to account low precision of hydrolakes
        hydrolake = gpd.overlay(hydrolake,extentpoly,how='intersection')        # Clip to extent in correct EPSG, ensure clipped area is correct
    hydrolake.to_file("%s%s_HydroLAKES.shp" %(folders[0],AOI))
    watermaskname = 'none'
