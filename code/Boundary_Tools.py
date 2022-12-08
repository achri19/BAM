#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 24, 2022
@author: Alchrist

"""
##############################################################################
import rasterio
import numpy as np
import os
from shapely.geometry import Polygon, LineString
#from shapely.ops import cascaded_union
#from shapely.ops import polygonize, unary_union
#from shapely.validation import make_valid

import geopandas as gpd
import warnings; warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)
import pandas as pd
pd.options.mode.chained_assignment = None
from rasterstats import zonal_stats

import warnings

#from terracatalogueclient import Catalogue
from tqdm.auto import tqdm  # provides a progressbar

import sys
if os.getcwd().split('/')[1] == 'Users':
    code_path = '/users/alchrist/documents/github/ANUGA/processing/code/'
else:
    code_path = '/scratch_lg/loac_hydro/alchrist/code'
sys.path.insert(1, code_path)
from make_tides import maketides

############################################################################
############################################################################
############################################################################
## PART VIII
## OUTPUT: Boundary condition information
##########################################################################

def set_boundary_conditions(delta,folders,res,parameters):
    step = 8
    print('\n\n\n##############################################################################################')
    print('##############################[Step %s][Set_Boundary_Conditions]###############################'%(step))
    print('##############################################################################################\n')
    ###########################################################################
    ###################### Import Config Parameters ###########################
    ###########################################################################
    xres,yres = res,res
    EPSG = parameters['EPSG'][0]                                                # Coordinate System must be UTM

    extentpoly = gpd.read_file("%s%s_modeldomain.shp" %(folders[7],delta))
    ulx,lry,lrx,uly = extentpoly.total_bounds   # Coordinates converted to UTM coordinate system
    extentpoly_line = extentpoly.geometry.boundary

    ###############################################################################
    ######################### Set final files for model ###########################
    # Set Boundaries. Open, tidal boundary determined as boundary with lowest elevation.
    # Other boundaries are set as transmissive Bi2
    # Set Upstream Boundary Condition at largest river

    ## Boundaries: North = 0, East = 1, South = 2, West = 3
    print('\n[Step %s][Set_Boundary_Conditions][Determine boundary type of each model side] .......\n'%(step))
    # sides = {'id':[0,1,2,3],'geometry':[LineString([(ulx,uly),(lrx,uly)]),LineString([(lrx,uly),(lrx,lry)]),LineString([(ulx,lry),(lrx,lry)]),LineString([(ulx,uly),(ulx,lry)])]}
    # df_line = gpd.GeoDataFrame(sides,columns=['id','geometry'])
    # df_line['geometry'] = df_line.geometry.buffer(1000)
    # df_line.crs = extentpoly.crs
    # boundaries = ['BoundaryType1','BoundaryType1','BoundaryType1','BoundaryType1'] # Boundary types: Bi tides, Bi2 transmissive, but no stage set

    x,y = extentpoly_line.geometry[0].coords.xy
    xy = pd.DataFrame(list(zip(x,y)), columns=['LON', 'LAT'])

    geoms = []
    ids = []
    for xy1 in range(len(xy)-1):
        x1 = xy.iloc[xy1]['LON']
        y1 = xy.iloc[xy1]['LAT']
        x2 = xy.iloc[xy1+1]['LON']
        y2 = xy.iloc[xy1+1]['LAT']
        geoms.append(LineString([(x1,y1),(x2,y2)]))
        ids.append(xy1)
    sides = {'id':ids,'geometry':geoms}
    df_line = gpd.GeoDataFrame(sides,columns=['id','geometry'])
    df_line['geometry'] = df_line.geometry.buffer(500)
    df_line.crs = extentpoly.crs
    df_line['index'] = df_line.index

    print('\n[Step %s][Set_Boundary_Conditions][Largest river channel along model edge = upstream boundary, river discharge conditions] .......\n'%(step))
    rivers = gpd.read_file("%s%s_rivers_%s.shp" %(folders[7],delta,xres))
    rivers = gpd.overlay(rivers,extentpoly,how='intersection')

    intersections = gpd.overlay(df_line,rivers,how='intersection')
    intersections['mean'] = pd.DataFrame(zonal_stats(vectors = intersections['geometry'],raster = "%s/%s_%s_%s.tif" %(folders[4] + elev_name,delta,elev_name,xres),stats='mean',nodata = '-9999'))['mean']
    intersections['areadepth'] = intersections.area * intersections['mean']
    try:
        #inlet = intersections.loc[intersections['mean']==min(intersections['mean'])].iloc[0]
        #inlet = intersections.loc[intersections.area==max(intersections.area)]
        inlet = intersections.loc[intersections['areadepth'] == min(intersections['areadepth'])]
    except:
        print('no intersections between river and model boundary')
        upstreamX = 0
        upstreamY = 0
    else:
        upstreamX= int(inlet['geometry'].centroid.x)
        upstreamY= int(inlet['geometry'].centroid.y)
        print('##################### Discharge boundary conditions set at %s,%s' %(upstreamX,upstreamY))

    #boundaries = ['BoundaryType1'] * (len(df_line)-1) # Boundary types: Bi tides, Bi2 transmissive, but no stage set
    df_line['boundary'] = ['BoundaryType1'] * (len(df_line))
    if os.path.isfile('%s%s_tidebnd.shp' %(folders[0],delta)) == True:
        print('\n[Step %s][Set_Boundary_Conditions][%s_tidebnd.shp will determine downstream boundary, tidal conditions] .......\n'%(step,delta))
        tide_bnd = gpd.read_file('%s%s_tidebnd.shp' %(folders[0],delta))
        tide_bnd = tide_bnd.to_crs(df_line.crs)
        tideboundary = gpd.overlay(df_line,tide_bnd,how='intersection')
        tide = tideboundary['index']
        #b#oundaries = np.where(boundaries[.index].isin(tide),'BoundaryType2',boundaries)
        #boundaries = np.where(df_line.iloc[tide],'BoundaryType2','BoundaryType1')
        df_line['boundary'].loc[tide] = 'BoundaryType2'

        tide_bnd = tide_bnd['geometry'].centroid
        tide_bnd = tide_bnd.to_crs("EPSG:4326")
        tidex= float(tide_bnd.x)
        tidey= float(tide_bnd.y)
        if tidex < 0:
            tidex = 360+tidex
    else:
        print('\n[Step %s][Set_Boundary_Conditions][Deepest model side = downstream boundary, tidal conditions] .......\n'%(step))
        df_line['mean'] = pd.DataFrame(zonal_stats(vectors = df_line['geometry'],raster = "%s%s_bathy_%s.tif" %(folders[8],delta,xres),stats='mean',nodata = '-9999'))['mean']
        print('\n[Step %s][Set_Boundary_Conditions][Deepest model side = downstream boundary, tidal conditions] .......\n'%(step))
        tide = df_line.loc[df_line['mean']==min(df_line['mean'])]
        df_line['boundary'].loc[tide] = 'BoundaryType2'
        #deepest = int(df_line['id'][df_line['mean']==min(df_line['mean'])])

        ##############################################################################
        ##############################################################################
        oceans = gpd.read_file("%s%s_fulloceans_%s.shp" %(folders[7],delta,xres))
        oceans = oceans.to_crs(df_line.crs)
        oceanboundary = gpd.overlay(df_line,oceans,how='intersection')
        oceanboundary['mean'] = pd.DataFrame(zonal_stats(vectors = oceanboundary['geometry'],raster = "%s%s_bathy_%s.tif" %(folders[8],delta,xres),stats='mean',nodata = '-9999'))['mean']
        tideboundary = oceanboundary.loc[oceanboundary['mean']==min(oceanboundary['mean'])]
        tideboundary = tideboundary['geometry'].centroid
        tideboundary = tideboundary.to_crs("EPSG:4326")
        tidex= float(tideboundary.x)
        tidey= float(tideboundary.y)
        if tidex < 0:
            tidex = 360+tidex
    ## Get Tide Data
    ## Run pyTMD to get global tidal predictions
    ## Set downstream boundary conditions
    ##############################################################################
    ##############################################################################
    ##############################################################################
    print('\n[Step %s][Set_Boundary_Conditions][Set model run period] .......\n'%(step))
    startdate = '20100101'
    enddate = '20211001'
    print('##################### Default simulation start and end are %s - %s' %(startdate, enddate))
    print('\n[Step %s][Set_Boundary_Conditions][Get water stage time series from TPXO Global Tide Model] .......\n'%(step))

    try:
        boundary_data = np.genfromtxt("%s/%s_tides_lat%s_lon%s_%s.csv" %(folders[5],delta,np.round(tidey,2),np.round(tidex,2),startdate),delimiter=',')
    except:
        boundary_data = np.column_stack((maketides(float(tidey),float(tidex),str(startdate),str(enddate),10)))# value every 10 minutes
        np.savetxt("%s/%s_tides_lat%s_lon%s_%s.csv" %(folders[5],delta,np.round(tidey,2),np.round(tidex,2),startdate),boundary_data,delimiter=',')
    print('##################### Using tide data from file: %s_tides_lat%s_lon%s_%s.csv' %(delta,np.round(tidey,2),np.round(tidex,2),startdate))

    print('\n[Step %s][Set_Boundary_Conditions] Finished .......\n'%(step))

    return df_line, upstreamX,upstreamY,tidex,tidey



def get_inlet(delta,folders,res,parameters,inletx,inlety,elev_name):
    step = '8b'
    print('\n\n\n##############################################################################################')
    print('##############################[Step %s][Get Inlet Location]###############################'%(step))
    print('##############################################################################################\n')
    ###########################################################################
    ###################### Import Config Parameters ###########################
    ###########################################################################
    xres,yres = res,res
    EPSG = parameters['EPSG'][0]                                                # Coordinate System must be UTM

    extentpoly = gpd.read_file("%s%s_modeldomain.shp" %(folders[7],delta))
    ulx,lry,lrx,uly = extentpoly.total_bounds   # Coordinates converted to UTM coordinate system
    extentpoly_line = extentpoly.geometry.boundary
    print('\n[Step %s][Get Inlet Location][Find River Boundary] .......\n'%(step))


    x,y = extentpoly_line.geometry[0].coords.xy
    xy = pd.DataFrame(list(zip(x,y)), columns=['LON', 'LAT'])

    geoms = []
    ids = []
    for xy1 in range(len(xy)-1):
        x1 = xy.iloc[xy1]['LON']
        y1 = xy.iloc[xy1]['LAT']
        x2 = xy.iloc[xy1+1]['LON']
        y2 = xy.iloc[xy1+1]['LAT']
        geoms.append(LineString([(x1,y1),(x2,y2)]))
        ids.append(xy1)
    sides = {'id':ids,'geometry':geoms}
    df_line = gpd.GeoDataFrame(sides,columns=['id','geometry'])
    df_line['geometry'] = df_line.geometry.buffer(500)
    df_line.crs = extentpoly.crs
    df_line['index'] = df_line.index
    if (inlety == -9999) & (inletx == -9999) :

        print('\n[Step %s][Get Inlet Location][Largest river channel along model edge = upstream boundary, river discharge conditions] .......\n'%(step))
        rivers = gpd.read_file("%s%s_rivers_%s.shp" %(folders[7],delta,xres))
        rivers = gpd.overlay(rivers,extentpoly,how='intersection')

        intersections = gpd.overlay(df_line,rivers,how='intersection')
        intersections['mean'] = pd.DataFrame(zonal_stats(vectors = intersections['geometry'],raster = "%s/%s.tif" %(folders[4],elev_name),stats='mean',nodata = '-9999'))['mean']
        intersections['areadepth'] = intersections.area * intersections['mean']
        intersections.to_file('test2.shp')
        try:
            #inlet = intersections.loc[intersections['mean']==min(intersections['mean'])].iloc[0]
            #inlet = intersections.loc[intersections.area==max(intersections.area)]
            inlet = intersections.loc[intersections['areadepth'] == np.nanmin(intersections['areadepth'])]
        except:
            print('no intersections between river and model boundary')
            inletx = 0
            inlety = 0
        else:
            inletx= int(inlet['geometry'].centroid.x)
            inlety= int(inlet['geometry'].centroid.y)
            print('##################### Discharge boundary conditions set at %s,%s' %(inletx,inlety))

    #boundaries = ['BoundaryType1'] * (len(df_line)-1) # Boundary types: Bi tides, Bi2 transmissive, but no stage set

    print('\n[Step %s][Get Inlet Location] Finished .......\n'%(step))

    return df_line, round(inletx,0),round(inlety,0)


def get_tidal_boundary(delta,folders,res,parameters,tide_bnd,tidey,tidex):
    step = '8A'
    print('\n\n\n##############################################################################################')
    print('##############################[Step %s][Get Tidal Boundary]###############################'%(step))
    print('##############################################################################################\n')
    ###########################################################################
    ###################### Import Config Parameters ###########################
    ###########################################################################
    xres,yres = res,res
    EPSG = parameters['EPSG'][0]                                                # Coordinate System must be UTM

    extentpoly = gpd.read_file("%s%s_modeldomain.shp" %(folders[7],delta))
    ulx,lry,lrx,uly = extentpoly.total_bounds   # Coordinates converted to UTM coordinate system
    extentpoly_line = extentpoly.geometry.boundary

    print('\n[Step %s][Get Tidal Boundary][Find Tidal Boundary] .......\n'%(step))
    x,y = extentpoly_line.geometry[0].coords.xy
    xy = pd.DataFrame(list(zip(x,y)), columns=['LON', 'LAT'])

    geoms = []
    ids = []
    for xy1 in range(len(xy)-1):
        x1 = xy.iloc[xy1]['LON']
        y1 = xy.iloc[xy1]['LAT']
        x2 = xy.iloc[xy1+1]['LON']
        y2 = xy.iloc[xy1+1]['LAT']
        geoms.append(LineString([(x1,y1),(x2,y2)]))
        ids.append(xy1)
    sides = {'id':ids,'geometry':geoms}
    df_line = gpd.GeoDataFrame(sides,columns=['id','geometry'])
    df_line['geometry'] = df_line.geometry.buffer(500)
    df_line.crs = extentpoly.crs
    df_line['index'] = df_line.index

    df_line['boundary'] = ['NonTidal'] * (len(df_line))
    df_line.to_file('test.shp')
    if len(tide_bnd) > 0:
        df_line_tide = gpd.sjoin(df_line,tide_bnd, how='inner', op='within')
        df_line_tide['id_left'].values
        df_line['boundary'][df_line_tide['id_left'].values] = 'Tide'
    else:
        df_line['mean'] = pd.DataFrame(zonal_stats(vectors = df_line['geometry'],raster = "%s%s_bathy_%s.tif" %(folders[8],delta,xres),stats='mean',nodata = '-9999'))['mean']
        print('\n[Step %s][Set_Boundary_Conditions][Deepest model side = downstream boundary, tidal conditions] .......\n'%(step))
        # tide = df_line.loc[df_line['mean']==min(df_line['mean'])]
        df_line['boundary'][df_line['mean']==np.nanmin(df_line['mean'])] = 'Tide'
        #deepest = int(df_line['id'][df_line['mean']==min(df_line['mean'])])
    if (tidex == -9999) & (tidey == -9999) :
        tide_lines = df_line[df_line['boundary']=='Tide'].geometry
        #tide_lines = tide_lines.to_crs('EPSG:4326')

        tidex = tide_lines[tide_lines.length == np.nanmax(tide_lines.length)].centroid.iloc[0].coords.xy[0][0]
        tidey = tide_lines[tide_lines.length == np.nanmax(tide_lines.length)].centroid.iloc[0].coords.xy[1][0]
    print('\n[Step %s][Get Tidal Boundary] Finished .......\n'%(step))

    return df_line,round(tidex,2),round(tidey,2)


def get_tide_data_pytmd(delta,path,tidex,tidey):
    step = '9'
    if tidex < 0:
        tidex = 360+tidex
    ## Get Tide Data
    ## Run pyTMD to get global tidal predictions
    ## Set downstream boundary conditions
    ##############################################################################
    ##############################################################################
    ##############################################################################
    print('\n[Step %s][Set_Boundary_Conditions][Set model run period] .......\n'%(step))
    startdate = '20100101'
    enddate = '20211001'
    print('##################### Default simulation start and end are %s - %s' %(startdate, enddate))
    print('\n[Step %s][Set_Boundary_Conditions][Get water stage time series from TPXO Global Tide Model] .......\n'%(step))

    try:
        tide_data = np.genfromtxt("%s/%s_tides_lat%s_lon%s_%s.csv" %(path,delta,np.round(tidey,2),np.round(tidex,2),startdate),delimiter=',')
    except:
        tide_data = np.column_stack((maketides(float(tidey),float(tidex),str(startdate),str(enddate),10)))# value every 10 minutes
        np.savetxt("%s/%s_tides_lat%s_lon%s_%s.csv" %(path,delta,np.round(tidey,2),np.round(tidex,2),startdate),tide_data,delimiter=',')
    print('##################### Using tide data from file: %s_tides_lat%s_lon%s_%s.csv' %(delta,np.round(tidey,2),np.round(tidex,2),startdate))

    print('\n[Step %s][Set_Boundary_Conditions] Finished .......\n'%(step))

    return tide_data
