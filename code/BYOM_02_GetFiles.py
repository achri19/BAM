#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 07:28:56 2021

@author: Alchrist
Build Your Own Model
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
    path = '/Volumes/FortressL3/ANUGA/GlobalDeltas/'
    path_code = '/Users/Alchrist/Documents/Github/ANUGA/processing/code/'
    #build_path = '/Users/Alchrist/Documents/Github/globalcoast/'
    path_templates = '/Users/Alchrist/Documents/Github/globalcoast/templates/'
    path_configs = '/Users/Alchrist/Documents/Github/globalcoast/configs/'
    cnes_tools_path = '/Volumes/FortressL3/Temp_ANUGA/Tools/swot-hydroloy-toolbox/'
    path_ancillary = path + 'inputs/'
    path_examples = path + 'examples/'
else:
    path = '/projects/loac_hydro/alchrist/anuga/'
    path_code = path + 'code/'
    #build_path = '/scratch_lg/loac_hydro/alchrist/anuga/build/'
    path_templates = path + 'templates/'
    path_configs = path + 'configs/'
    path_ancillary = path + 'ancillary/'
    cnes_tools_path = '/scratch_sm/loac_hydro/swot-hydrology-toolbox/'
    #mangroves = gpd.read_file(path + "inputs/Temp_ANUGA/GMW_2016_v2_fixed.shp")
    path_examples = path + '/examples/'

sys.path.insert(1,path_code)
print(path_code)
from BYOM_Utilities_V1 import (build_directory,
                               setup_AOI_files,
                               make_polygons,
                               make_channel_networks,
                               make_model_foundation,
                               set_boundary_conditions,
                               make_watermask)

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
res = 30

configfile = path_configs + 'Deltas_SetupParameters_MASTER_V2.csv'
parameter_file = pd.read_csv(configfile)
AOIs = parameter_file['AOI']

i=0

for AOI in AOIs[1:]:
    skip = False
    print('################ ################# ##################')
    print('################ ################# ##################')
    print('################ ################# ##################')
    print('################ ANUGA Model Setup ##################')
    print('################ ################# ##################')
    print('################ ################# ##################')
    print('################ ################# ##################')

    print('Before building the model, we need 3 pieces: \n1) model domain (shapefile) \n2) approximate tide boundary (shapefile) \n3) watermask')

    AOI = AOI.lower()
    parameters = parameter_file[parameter_file['AOI']==AOI.capitalize()].reset_index(drop=True)

    ## STEP 1
    ## Set up working directory and build sub folders if needed
    print('\n')
    print('Study area is ' + AOI)
    print('Resolution of this setup is %sm' %(res))
    working_path,folders = build_directory(path_examples, AOI)

    final_watermasks = [os.path.join(dirpath,f)
        for dirpath,dirnames, files in os.walk(folders[0])
        for f in fnmatch.filter(files,'*_finalwatermask*.tif')]

    if len(final_watermasks) >0:
        print('You already have a chosen water mask:')
        print(final_watermasks)

        ref_10m,parameters = get_extent_parameters(path_ancillary,AOI,folders,res,parameters)

        ref = setup_AOI_files(path,
                            AOI,
                            folders,
                            res,
                            parameters)


        parameters = pd.read_csv('%s/%s_Configuration.csv' %(folders[2],AOI))
        print(parameters.iloc[0] )

        EPSG = parameters['EPSG'][0]
        ulx = parameters['ulx'][0]
        uly = parameters['uly'][0]
        lrx = parameters['lrx'][0]
        lry = parameters['lry'][0]

        model_domain = gpd.read_file('%s%s_modeldomain.shp' %(folders[7],AOI))
        AOI_extent = gpd.read_file('%s%s_extent_%s.shp' %(folders[7],AOI,EPSG))

        clean_with_landcover = False
        ref = rasterio.open('%s/%s_GEBCO_%s.tif' %(folders[8],AOI,res))
        watermaskname = make_watermask(path_ancillary,
                                       AOI,
                                       folders,
                                       parameters,
                                       ref,clean_with_landcover,
                                       os.path.isfile('%s%s_landmask_%s.tif' %(folders[8],AOI,res)))

        watermask = rasterio.open('%s%s_watermask_%s.tif' %(folders[8],AOI,res)).read(1)

        make_polygons(AOI,
                        folders,
                        parameters,
                        ref,
                        watermaskname,
                        path_templates,
                        os.path.isfile("%s%s_farocean_%s.shp" %(folders[7],AOI,res)))


        segment_width = 30
        pixel_step = int(round(segment_width/res))
        segment_width_2 = 60
        pixel_step_2 = int(round(segment_width_2/res))
        try:
            distance,widths = make_channel_networks(folders,
                                                  AOI,
                                                  ref,
                                                  parameters,
                                                  pixel_step,
                                                  (os.path.isfile("%s%s_widths_%sx%s.tif" %(folders[8],AOI,res,pixel_step))) | (os.path.isfile("%s%s_widths_%sx%s.tif" %(folders[8],AOI,res,pixel_step_2))))
        except:
            distance,widths = make_channel_networks(folders,
                                                  AOI,
                                                  ref,
                                                  parameters,
                                                  pixel_step_2,
                                                  os.path.isfile("%s%s_widths_%sx%s.tif" %(folders[8],AOI,res,pixel_step_2)))
            pixel_step = int(round(segment_width_2/res))

        elevation,elev_name,elevationpath = make_model_foundation(path,
                                                        parameters,
                                                        AOI,
                                                        folders,
                                                        ref,
                                                        distance,
                                                        widths,
                                                        watermask,pixel_step,path
                                                        )

        ## STEP 8
        boundaries,upstreamX,upstreamY,tideX,tideY = set_boundary_conditions(AOI,folders,res,parameters)

        ## STEP 5
        ## Make river and ocean segments for SWOT simulator
        ## Approximately 900m for river segments and 10000m2 for oceans
        # swot_res = res*3
        # try:
        # dem = rasterio.open('%s/%s_%s.tif' %(elevationpath, AOI, elev_name))
        # superpixel_area = 10000
        # #ref = rasterio.open('%s%s_%s.tif' %(dempath,AOI,dempath.split('/')[-2]))
        # make_segments_for_swot(folders,AOI,ref,parameters,superpixel_area,False)
        #

        #
        # #############################################################################
        # #############################################################################
        # # Make the script to build the model (must be run in Python ANUGA env)
        # #############################################################################
        # #############################################################################
        #
        # now = datetime.now()
        #
        # filein = open(build_path + '/templates/Build_HPCModel_V1.py')
        # replacements = {'res':res,
        #                 'config': 'config_%s.csv' %(AOI),#,datetime.now().strftime("%Y%m%d")),
        #                 'SWOT_HYDROLOGY_TOOLBOX':'$SWOT_HYDROLOGY_TOOLBOX'}
        # template = Template(filein.read())
        # print('\n#####################  Saving run script as %s_Build_HPCModel_V1.py\n' %(AOI))
        # runanuga = template.substitute(replacements)
        # file = open('%s/%s_Build_HPCModel_V1.py' %(elevationpath,AOI),'w')
        # file.write(runanuga)
        # file.close()
        #
        # ##############################################################################
        # ##############################################################################
        # ##############################################################################
        # discharges = pd.read_csv(build_path + '/configs/Deltas_Discharge.csv')
        # runanuga = pd.DataFrame(columns = ['AOI','Path','HPCPath','EPSG','TideX','TideY','BaseScale','UpstreamX','UpstreamY','BoundaryType1','BoundaryType2','North','East','South','West','Discharge_cms','ElevationName','ManningLUT','res'])
        # runanuga.loc[0,'AOI'] = delta
        # #runanuga.loc[0,'Path'] = deltapath
        # #runanuga.loc[0,'HPCPath'] = '/scratch_lg/loac_hydro/alchrist/anuga/build/examples/%s/' %(delta)
        # runanuga.loc[0,'TideX'] = tideX
        # runanuga.loc[0,'TideY'] = tideY
        # runanuga.loc[0,'ElevationName'] = elev_name
        # runanuga.loc[0,'Discharge_cms'] = discharges[discharges['AOI']==delta.capitalize()].reset_index(drop=True)['HydroSHEDS_Avg_Discharge'][0].astype('int')
        # runanuga.loc[0,'EPSG'] = parameters['EPSG'][0]
        # runanuga.loc[0,'LandcoverMethod'] = parameters['LandcoverMethod'][0]
        # runanuga.loc[0,'BoundaryType1'] = 'Br'
        # runanuga.loc[0,'BoundaryType2'] = 'Bout'
        # runanuga.loc[0,'North'] = boundaries[0]
        # runanuga.loc[0,'East']  = boundaries[1]
        # runanuga.loc[0,'South'] = boundaries[2]
        # runanuga.loc[0,'West']  = boundaries[3]
        # runanuga.loc[0,'AOI'] = delta
        # runanuga.loc[0,'UpstreamX'] = round(upstreamX)
        # runanuga.loc[0,'UpstreamY'] = round(upstreamY)
        # runanuga.loc[0,'StartDate'] = '20210101'
        # runanuga.loc[0,'EndDate'] = '20210107'
        # # runanuga['Mesh'] = 'base'
        # # runanuga['Scale'] = 10000
        # runanuga['Mesh'] = 'base'
        # runanuga['Scale'] = min(riverscale,10000)
        #
        # runanuga['Mesh'] = 'river'
        # runanuga.loc[1,'Mesh'] =  'fullocean'
        # runanuga.loc[2,'Mesh'] =  'land'
        # runanuga.loc[3,'Mesh'] =  'lake'
        # runanuga['Scale'] = min(riverscale,10000)
        # runanuga.loc[1,'Scale'] =  100000
        # runanuga.loc[2,'Scale'] =  20000
        # runanuga.loc[3,'Scale'] =  50000
        #
        # runanuga.loc[0,'ulx']  = parameters['ulx'][0]                                                 # ULX coordinate
        # runanuga.loc[0,'lry']  = parameters['lry'][0]                                                 # LRY coordinate
        # runanuga.loc[0,'lrx']  = parameters['lrx'][0]                                                 # LRX coordinate
        # runanuga.loc[0,'uly']  = parameters['uly'][0]
        # runanuga.loc[0,'res']  = res
        #
        # runanuga.to_csv('%s/config_%s.csv' %(elevationpath,delta))#,datetime.now().strftime("%Y%m%d")))
        # print('To adjust mesh type, boundary types, and mesh resolution, make changes to %s/config_%s.csv' %(folders[4],delta))#,datetime.now().strftime("%Y%m%d")))

        print('Cleaning up temporary files')
        try:shutil.rmtree(folders[1])
        except:''
        i = i+1
