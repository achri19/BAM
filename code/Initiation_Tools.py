
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 24, 2022
@author: Alchrist

Initiate Tools includes functions for choosing AOI, building a new directory, making extent file, getting EPSG and boundary locations.

"""

import fnmatch
import os
import numpy as np
import rasterio
from shapely.geometry import Polygon, LineString
import geopandas as gpd
import math

def build_directory(input_path, aoi):
    '''
    Parameters
    ----------
    input_path : string
        root directory of files (inputs and GEE watermasks).
    aoi : string
        AOI.
    '''
    print('\n\n\n##############################################################################################')
    print('################################[Step 1][Build Directory]#####################################')
    print('##############################################################################################\n')

    try:
        aoipath = [os.path.join(dirpath,f)
                for dirpath,dirnames, files in os.walk(input_path)
                for f in fnmatch.filter(dirnames,'%s' %(aoi))][0] +'/'
    except:
        aoipath = input('The directory %s does not exist, what is the working directory:')
    else:
        print('##################### The working directory set as: \n\n%s\n ' %(aoipath))

    ## Build model directory
    tier1_folders = [aoipath + x for x in ['User_Defined_Files/','tmp/','Setup_Files/','Meshes/','DEMs/','Boundaries','Simulations/']]
    setup_path = tier1_folders[2]
    tier2_folders = [setup_path + x for x in['Setup_SHP/','Setup_RST/','Setup_FIG/']]
    folders = tier1_folders + tier2_folders
    print('##################### Folders are:')
    print('##################### 0 User_Defined_Files --> User shapefile of model domain and water mask')
    print('##################### 1 tmp --> For temporary files')
    print('##################### 3 Meshes --> Where we will build model meshes')
    print('##################### 4 DEMs --> Where we will build digital elevation models')
    print('##################### 5 Boundaries --> Where we will store boundary files')
    print('##################### 6 Simulations --> Where we will run simulations')
    print('##################### 7 Setup_Files/Setup_SHP --> Shapefiles for setup ')
    print('##################### 8 Setup_Files/Setup_RST --> Rasters for setup')
    print('##################### 9 Setup_Files/Setup_FIG --> Figures')

    for folder in folders:
        try: os.mkdir(folder)
        except: ''

    print('\n[Step 1][Build Directory] Finished .......\n')
    return aoipath,folders

############################################################################
############################################################################
############################################################################
############################################################################
## To determine the boundary coordinates for the AOI and get the correction
## EPSG coordinate system. Must be in UTM for ANUGA models
## INPUT: Delta name, working directory, folders
## OUTPUT: extent, bathymetry, topography, hydropolys, landcover, ndwi watermask
##########################################################################
def get_extent_parameters(input_path,aoi,folders,xres,parameters):
    '''
    Parameters
    ----------
    path : string
        root directory of files (inputs and GEE watermasks).
    aoi : string
        AOI name, must match input shapefile and folder.
    folders : np.array of strings
        folders within deltapath for model files.
    parameters : pd.dataframe
        Configuration parameters for model setup.
    '''
    step = 2
    print('\n\n\n##############################################################################################')
    print('#################################[Step %s][Setup_AOI_Files]####################################' %(step))
    print('##############################################################################################\n')

    # Set model resolution, set methods for building DEM of land and ocean, and set method for land cover classification
    ref_res = 10

    ##############################################################################
    ############################## Get AOI Extent ################################
    ##############################################################################

    # Model extent is defined with shapefile with naming format AOI_input.shp
    # Must be stored in User_Defined_Files folder within working directory
    print('\n[Step %s][Setup_AOI_Files][AOI extent] .......\n' %(step))
    print('##################### AOI extent is set by: %s_input.shp' %(aoi))

    if os.path.isfile('%s%s_input.shp' %(folders[0],aoi)):
        AOI = gpd.read_file('%s%s_input.shp' %(folders[0],aoi))
    else:
        inputfile = input('Model extent file %s%s_input.shp does not exist. What is the full path of model extent file?' %(folders[0], aoi))
        AOI = gpd.read_file(inputfile)

    if AOI.crs is None:
        AOI.crs = 'EPSG:4326'
    if AOI.crs != 'EPSG:4326':
        print("The input shapefile is not in the correction projection (EPSG 4326), reprojecting to EPSG 4326")
        AOI = AOI.to_crs(4326)

    ##############################################################################
    ################ Get projected coordinate system for AOI #####################
    ##############################################################################

    ## Total bounding coordinates in WGS84 coordinates
    ulx_wgs,lry_wgs,lrx_wgs,uly_wgs = AOI.total_bounds

    ## Determine the UTM coordinate system
    ## ANUGA models are assumed in UTM
    ## identify north/south and east/west coordiantes
    x1,y1,x2,y2 = math.floor(ulx_wgs),math.floor(uly_wgs),math.floor(lrx_wgs),math.floor(lry_wgs)

    print('\n[Step %s][Setup_AOI_Files][Determine EPSG code and UTM zone] .......\n'%(step))
    print('##################### ANUGA Models must be in UTM')
    zone = int(np.ceil((ulx_wgs + 180)/6))

    if y1>=0 and y2>=0:
        NS = 'n'
        EPSG = 32200+zone
    elif y1>=0 and y2<0:
        NS = 'n'
        NS2 = 's'
        EPSG = 32200+zone
    else:
        NS = 's'
        y = abs(y1)
        EPSG = 32700+zone
    if x1>=0:
        EW = 'e'
    elif x1<0 and x2>=0:
        EW = 'e'
        EW2 = 'w'
    else:
        EW = 'w'
        x = abs(x1)
    print('##################### UTM Zone: %s%s' %(zone,NS))

    print('##################### EPSG: %s ' %(EPSG))
    AOI_warped = AOI.to_crs('EPSG:%s' %(EPSG))
    AOI_warped.to_file('%s%s_modeldomain.shp' %(folders[7],aoi))

    ## Get bounding coordinates in UTM
    xs1,ys1,xs2,ys2 = AOI_warped.total_bounds

    ## Determine north/south and east/west
    print('\n[Step %s][Setup_AOI_Files][Extending AOI by 1000m] .......\n' %(step))
    ulx = int(min(xs1,xs2) - 1000)
    uly = int(max(ys1,ys2) + 1000)
    lrx = int(max(xs1,xs2) + 1000)
    lry = int(min(ys1,ys2) - 1000)
    if (lrx-ulx)%xres !=0:
        lrx = int(lrx - (lrx-ulx)%xres)
    if (uly-lry)%xres !=0:
        uly = int(uly - (uly-lry)%xres)

    # Model domain extent
    print('\n[Step %s][Setup_AOI_Files][Setting up AOI extent] .......\n' %(step))
    extent = [(ulx),(uly)],[(lrx),(uly)],[(lrx),(lry)],[(ulx),(lry)]
    poly = Polygon([[p[0], p[1]] for p in extent])
    boundingbox = gpd.GeoDataFrame(index=[0],crs='EPSG:%s'%(EPSG),geometry=[poly])

    np.savetxt('%s%s_extent.csv' %(folders[0],aoi), extent,delimiter=',', fmt= '%1.3f') ## USED BY ANUGA MODEL
    boundingbox.to_file("%s%s_extent_%s.shp" %(folders[7],aoi,EPSG))
    extentarea = np.round(boundingbox.area[0],-6)
    print('##################### AOI bounds are : %s, %s, %s, %s' %(round(ulx,-1),round(lry,-1),round(lrx,-1),round(uly,-1)))
    print('##################### Approximate area of AOI extent is %s km^2' %(extentarea/1000000))

    ##############################################################################
    ############# GEBCO (Bathymetry, relative to mean sea level)  ################
    ##############################################################################
    # GEBCO file will be reference for projection, origin, resolution
    # Get GEBCO data within the model domain
    print('\n[Step %s][Setup_AOI_Files][Downloading GEBCO Dataset as reference for projection, resolution, etc] .......\n' %(step))
    try:
        ref_10m = rasterio.open('%s%s_GEBCO_%s.tif' %(folders[8],aoi,ref_res))
    except:
        os.system('gdalwarp -overwrite -tr %s %s %sgebco_2020_geotiff/gebco_all.vrt %s%s_GEBCO_%s.tif '\
                  '-t_srs EPSG:%s -te %s %s %s %s -srcnodata -9999 -dstnodata -9999 -co COMPRESS=DEFLATE -q'
                  %(ref_res,ref_res,input_path,folders[8],aoi,ref_res,EPSG,ulx,lry,lrx,uly))
        ref_10m = rasterio.open('%s%s_GEBCO_%s.tif' %(folders[8],aoi,ref_res))
    save_profile_10m = ref_10m.profile
    parameters['ulx'] = ulx
    parameters['lry'] = lry
    parameters['lrx'] = lrx
    parameters['uly'] = uly
    parameters['EPSG'] = EPSG
    parameters.to_csv('%s/%s_Configuration.csv' %(folders[2],aoi))

    print('\n[Step 2][Setup_AOI_Files][Saving configuration file] .......\n')
    print('##################### Saved as %s%s_Configuration.csv' %(folders[2].split(aoi)[-1],aoi))
    return ref_10m,parameters
