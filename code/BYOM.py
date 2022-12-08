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

if os.getcwd().split('/')[1] == 'Users':
    path = '/Users/alchrist/Documents/GitHub/ANUGA/'
    path_code = path + 'processing/code/'
    path_templates = path + 'processing/templates/'
    path_configs = path  + 'processing/configs/'
    cnes_tools_path = path + 'swot-hydroloy-toolbox/'
    path_ancillary = path + 'ancillary/'
    path_examples = '/Volumes/FortressL3/ANUGA/GlobalDeltas/SWOT_Sim/' #''/Volumes/FortressL3/ANUGA/GlobalDeltas/examples/'
else:
    path = '/projects/loac_hydro/alchrist/anuga/'
    path_code = path + 'code/'
    path_templates = path + 'templates/'
    path_configs = path + 'configs/'
    path_ancillary = path + 'ancillary/'
    cnes_tools_path = '/scratch_sm/loac_hydro/swot-hydrology-toolbox/'
    path_examples = path + '/examples/'

sys.path.insert(1,path_code)
from BYOM_Utilities_V1 import build_directory, get_extent_parameters, setup_AOI_files, make_polygons,make_channel_networks,make_model_foundation, set_boundary_conditions, make_watermask
from BYOM_extra import make_mesh_polygons, make_SWORD_networks, make_segments_for_swot
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

configfile = path_configs + '/Deltas_SetupParameters_MASTER_V2.csv'
parameter_file = pd.read_csv(configfile)
deltas = parameter_file['AOI']

i=0

deltas = ['gabon']
for delta in deltas[:]:
    skip = False
   
    print('Before building the model, we need 3 pieces: \n1) model domain (shapefile) \n2) approximate tide boundary (shapefile) \n3) watermask')

    gee_ndwi_path = path + 'GEE_NDWI_watermask' #'/Volumes/GoogleDrive/My Drive/GEE_drive/'
    print('The watermask should be stored in %s' %(gee_ndwi_path))
    print('If it is not there yet, we will make one. But this will end the program.')
    print('You will need to download the file from google drive when it is made and move it to this folder')

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

    ## Set up working directory and build sub folders if needed
    print('\n')
    print('Study area is ' + AOI)
    print('Resolution of this setup is %sm' %(res))
    working_path,folders = build_directory(path_examples, AOI)
    ref_10m,parameters = get_extent_parameters(path_ancillary,AOI,folders,res,parameters)

    ## STEP 2
    ## Download ocean bathyemtry, land topography, landcover,  and water mask datasets
    ref = setup_AOI_files(path,
                        delta,
                        folders,
                        res,
                        parameters)

    parameters = pd.read_csv('%s/%s_Configuration.csv' %(folders[2],delta))
    print(parameters.iloc[0] )

    EPSG = parameters['EPSG'][0]
    # gee = [os.path.join(dirpath,f)
    #     for dirpath,dirnames, files in os.walk(gee_ndwi_path)
    #     for f in fnmatch.filter(files,'%s*' %(delta))]
    gee = [os.path.join(dirpath,f)
        for dirpath,dirnames, files in os.walk(folders[0])
        for f in fnmatch.filter(files,'%s*_finalwatermask*' %(delta))]
    if len(gee) == 0:
         print('No watermask file exists, you must create on in GEE and download it')
         os.system('python %s/run_gee_V2.py %s %s%s_modeldomain.shp %s' %(code_path,delta,folders[6],delta,folders))
         print("We've made several options for water masks, you must download them and pick the best one before continuing")
         print('This will take a while, depending on how big your AOI is')
         print('Please wait to proceed until they are finsihed. Then download them and choose the best option')
         print('Once you choose the best one, put in into the folder %s' %(gee_ndwi_path))
         #user_input = input('Did you put it there? If yes, we can proceed if you enter "y". If not, the remaining code will fail ')
         # gee = [os.path.join(dirpath,f)
         #     for dirpath,dirnames, files in os.walk(gee_ndwi_path)
         #     for f in fnmatch.filter(files,'%s*' %(delta))]
         # print('Stop it, download the file and move it now')
         print('\n\n\n\n')
         #continue

    print('\n')
    print('Study area is ' + delta)
    print('Resolution of this setup is %sm' %(res))
    print('\n')




    ## STEP 3A
    # Clean up water GEE watermask, improve connectivity
    watermaskname = make_watermask(path, delta,
                    folders,
                    parameters,
                    ref,
                    os.path.isfile('%s%s_watermask_%s.tif' %(folders[7],delta,res)))

    ## STEP 3B
    ## Make polygons for land, ocean, lakes, rivers,
    make_polygons(delta,
                    folders,
                    parameters,
                    ref,
                    watermaskname,
                    templates_path,
                    os.path.isfile("%s%s_farocean_%s.tif" %(folders[7],delta,res)))

    ## STEP 4
    ## Get distance, widths, segments, and channel networks with 150m reaches
    pixel_step = int(round(150/res))
    distance,widths,riverscale,med_width,cellsperwidth,watermask = make_channel_networks(folders,
                                                                                          delta,
                                                                                          ref,
                                                                                          parameters,
                                                                                          pixel_step,True)
                                                                                          os.path.isfile("%s%s_widths_%s.tif" %(folders[7],delta,round(pixel_step*res))))

    pixel_step = 4
    sword_res = 10
    sword_ref = rasterio.open('%s%s_GEBCO_%s.tif' %(folders[7],delta,sword_res))

    watermaskname_sword = make_watermask(path, delta,
                    folders,
                    parameters,
                    sword_ref,
                    os.path.isfile('%s%s_watermask_%s.tif' %(folders[7],delta,sword_res)))

    if (os.path.isfile("%s%s_river_centerline_%s.shp" %(folders[6],delta,int(sword_res*pixel_step)))==False) & (os.path.isfile("%s%s_river_centerline_%s.shp" %(folders[6],delta,int(sword_res*(pixel_step+1))))==False):
        try: distance = make_SWORD_networks(folders,sword_ref,delta,parameters,pixel_step,watermaskname_sword,False)
        except:
            distance = make_SWORD_networks(folders,sword_ref,delta,parameters,pixel_step+1,watermaskname_sword,False)

    ## STEP 6
  #   meshes = make_mesh_polygons(folders,
  #                               delta,
  #                               res,
  #                               parameters,
  #                               med_width,
  #                               cellsperwidth,False)
  # #                              os.path.isfile( "%s%s_land_mesh_%s.shp" %(folders[3],delta,res)))
    ## STEP 7
    pixel_step = int(round(150/res))

    # widths = rasterio.open('%s%s_widthsfromproximity_%s.tif' %(folders[7],delta,int(pixel_step*res))).read(1)
    elevation,elev_name,elevationpath = make_model_foundation(path,
                                                    parameters,
                                                    delta,
                                                    folders,
                                                    ref,
                                                    distance,
                                                    widths,
                                                    watermask,pixel_step,build_path,True
                                                    )


    ## STEP 8
    boundaries,upstreamX,upstreamY,tideX,tideY = set_boundary_conditions(delta,folders,res,parameters)

    ## STEP 5
    ## Make river and ocean segments for SWOT simulator
    ## Approximately 900m for river segments and 10000m2 for oceans
    swot_res = 100
    pixel_step = int(round(90/swot_res))
    target_segment_size = 10000
    scale = 1
    sigma = 0.95
    try:
        dem = rasterio.open(elevationpath + '/' + delta + '_' + elev_name + '_' + str(swot_res) + '.tif')
    except:
        den = gdal.Warp(elevationpath + '/' + delta + '_' + elev_name + '_' + str(swot_res) + '.tif', elevationpath + '/' + delta + '_' + elev_name + '.tif',dstNodata=-9999,srcNodata = -9999,xRes=res,yRes=res,creationOptions  = ["COMPRESS=DEFLATE"]).ReadAsArray()
        dem = rasterio.open(elevationpath + '/' + delta + '_' + elev_name + '_' + str(swot_res) + '.tif')

    make_segments_for_swot(folders,delta,dem,parameters,pixel_step,target_segment_size,scale,sigma,os.path.isfile('%s%s_oceanriver_bathy_superpixels_%sm2_%s_%s_felz_SWOTb.shp'%(folders[6],delta,target_segment_size,scale,sigma)))#,10,0.01)


    ##############################################################################
    ##############################################################################
    ## Make the script to build the model (must be run in Python ANUGA env)
    ##############################################################################
    ##############################################################################

    now = datetime.now()
    # filein = open(build_path + '/templates/Build_Model_V1.py')
    # template = Template(filein.read())
    # print('\n#####################  Saving run script as %s_Build_Model_V1.py\n' %(delta))

    # runanuga = template.substitute(replacements)
    # file = open('%s/%s_Build_Model_V1.py' %(elevationpath,delta),'w')
    # file.write(runanuga)
    # file.close()

    filein = open(build_path + '/templates/Build_HPCModel_V1.py')
    replacements = {'res':res,
                    'config': 'config_%s.csv' %(delta),#,datetime.now().strftime("%Y%m%d")),
                    'SWOT_HYDROLOGY_TOOLBOX':'$SWOT_HYDROLOGY_TOOLBOX'}
    template = Template(filein.read())
    print('\n#####################  Saving run script as %s_Build_HPCModel_V1.py\n' %(delta))
    runanuga = template.substitute(replacements)
    file = open('%s/%s_Build_HPCModel_V1.py' %(elevationpath,delta),'w')
    file.write(runanuga)
    file.close()

    ##############################################################################
    ##############################################################################
    ##############################################################################
    discharges = pd.read_csv(build_path + '/configs/Deltas_Discharge.csv')
    runanuga = pd.DataFrame(columns = ['AOI','Path','HPCPath','EPSG','TideX','TideY','BaseScale','UpstreamX','UpstreamY','BoundaryType1','BoundaryType2','North','East','South','West','Discharge_cms','ElevationName','ManningLUT','res'])
    runanuga.loc[0,'AOI'] = delta
    #runanuga.loc[0,'Path'] = deltapath
    #runanuga.loc[0,'HPCPath'] = '/scratch_lg/loac_hydro/alchrist/anuga/build/examples/%s/' %(delta)
    runanuga.loc[0,'TideX'] = tideX
    runanuga.loc[0,'TideY'] = tideY
    runanuga.loc[0,'ElevationName'] = elev_name
    runanuga.loc[0,'Discharge_cms'] = discharges[discharges['AOI']==delta.capitalize()].reset_index(drop=True)['HydroSHEDS_Avg_Discharge'][0].astype('int')
    runanuga.loc[0,'EPSG'] = parameters['EPSG'][0]
    runanuga.loc[0,'LandcoverMethod'] = parameters['LandcoverMethod'][0]
    runanuga.loc[0,'BoundaryType1'] = 'Br'
    runanuga.loc[0,'BoundaryType2'] = 'Bout'
    runanuga.loc[0,'North'] = boundaries[0]
    runanuga.loc[0,'East']  = boundaries[1]
    runanuga.loc[0,'South'] = boundaries[2]
    runanuga.loc[0,'West']  = boundaries[3]
    runanuga.loc[0,'AOI'] = delta
    runanuga.loc[0,'UpstreamX'] = round(upstreamX)
    runanuga.loc[0,'UpstreamY'] = round(upstreamY)
    runanuga.loc[0,'StartDate'] = '20210101'
    runanuga.loc[0,'EndDate'] = '20210107'
    # runanuga['Mesh'] = 'base'
    # runanuga['Scale'] = 10000
    runanuga['Mesh'] = 'base'
    runanuga['Scale'] = min(riverscale,10000)

    runanuga['Mesh'] = 'river'
    runanuga.loc[1,'Mesh'] =  'fullocean'
    runanuga.loc[2,'Mesh'] =  'land'
    runanuga.loc[3,'Mesh'] =  'lake'
    runanuga['Scale'] = min(riverscale,10000)
    runanuga.loc[1,'Scale'] =  100000
    runanuga.loc[2,'Scale'] =  20000
    runanuga.loc[3,'Scale'] =  50000

    runanuga.loc[0,'ulx']  = parameters['ulx'][0]                                                 # ULX coordinate
    runanuga.loc[0,'lry']  = parameters['lry'][0]                                                 # LRY coordinate
    runanuga.loc[0,'lrx']  = parameters['lrx'][0]                                                 # LRX coordinate
    runanuga.loc[0,'uly']  = parameters['uly'][0]
    runanuga.loc[0,'res']  = res

    runanuga.to_csv('%s/config_%s.csv' %(elevationpath,delta))#,datetime.now().strftime("%Y%m%d")))
    print('To adjust mesh type, boundary types, and mesh resolution, make changes to %s/config_%s.csv' %(folders[4],delta))#,datetime.now().strftime("%Y%m%d")))

    print('Cleaning up temporary files')
    try:shutil.rmtree(folders[1])
    except:''
    i = i+1
