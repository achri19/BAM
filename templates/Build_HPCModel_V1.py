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
Split into separate Build and Run code: 10-26-2021
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
# Standard packages
import os
import time
import sys
import fnmatch
import numpy as np
import scipy
import scipy.interpolate
import pandas as pd
from datetime import datetime
from string import Template
import geopandas as gpd
from shapely.geometry import Polygon, Point
import rasterio
from pathlib import Path

## ANUGA packages
import anuga
from anuga.utilities.plot_utils import Make_Geotif
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.utilities import animate
from anuga import myid, numprocs, finalize, barrier
from anuga.parallel.parallel_inlet_operator import Parallel_Inlet_operator

## Plotting modules
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000

#mpl.rcParams['animation.ffmpeg_path'] = r'/Users/alchrist/Downloads/ffmpeg'
## If you get an error with producing animation plots, uncomment the line above and add directory for ffmpeg (r'/PATH/TO/ffmpeg')
import matplotlib.animation as animation
writer = animation.FFMpegWriter(fps=60)

if os.getcwd().split('/')[1] == 'Volumes':
    path = '/Volumes/FortressL3/ANUGA/GlobalDeltas/'
    build_path = '/Users/Alchrist/Documents/Github/globalcoast/'
    code_path = '/Users/Alchrist/Documents/Github/globalcoast/code/'
    templates_path = build_path + 'templates/'
    config_path = build_path + 'configure_files/'
    cnes_tools_path = '/Volumes/FortressL3/ANUGA/GlobalDeltas/swot-hydrology-toolbox/'
    local_path = path + 'examples/'
    hpc_path = '/scratch_lg/loac_hydro/alchrist/anuga/build/examples/'
else:
    path = '/scratch_lg/loac_hydro/alchrist/anuga/build/'
    build_path = '/scratch_lg/loac_hydro/alchrist/anuga/build/'
    code_path = '/scratch_lg/loac_hydro/alchrist/code/'
    templates_path = build_path + 'templates/'
    config_path = build_path + 'configure_files/'
    cnes_tools_path = '/scratch_lg/loac_hydro/alchrist/swot-hydrology-toolbox/'
    local_path = path + 'examples/'
    hpc_path = '/scratch_lg/loac_hydro/alchrist/anuga/build/examples/'

## My packages
sys.path.insert(1, code_path)
from fix_polygons import getpolygonpoints, removenearbypoints

print('######################################################################')
print('##############################################')
print('################################# Running ANUGA Model')
print('##############################################')
print('######################################################################')

##############################################################################
##############################################################################
##############################################################################
## Initiate timer
now = datetime.now()
nowtime = now.strftime('%Y%m%d')
time00 = time.time()
print('################################# Date and Time: %s ' %(now))

res = '$res'

##############################################################################
##############################################################################
##############################################################################
## Model Setup
##############################################################################
##############################################################################
##############################################################################

##############################################################################
## Get model parameters from parameter file
## To run a new model, adjust parameter file
## Configuration file is made in Global_Deltas.py
## Change as needed

parameters = pd.read_csv('$config',dtype=str,delimiter=',')
print('$config' + ' is the configuration file used for model setup')

delta           = parameters['AOI'][0]  .lower()                     # Name of the AOI
print('')
print('################################# Area of Interest:                %s' %(delta))
print('')
deltapath = local_path + delta + '/'

folders         = [deltapath + 'User_Defined_Files/',
                   deltapath + 'tmp/',
                   deltapath + 'Setup_Files/',
                   deltapath + 'Meshes/',
                   deltapath + 'Scenarios/',
                   deltapath + 'Models/',
                   deltapath + 'Setup_Files/Setup_SHP/',
                   deltapath + 'Setup_Files/Setup_RST/',
                   deltapath + 'Setup_Files/Setup_FIG/']

EPSG            = int(parameters['EPSG'][0])                         # Coordinate system
elevation       = parameters['ElevationName'][0]                     # Name of the elevation file
elevationpath   = folders[4] + elevation

startdate       = int(parameters['StartDate'][0])                    # Date to start simulation
enddate         = int(parameters['EndDate'][0])                      # Date to start simulation
boundarytype1   = parameters['BoundaryType1'][0]                     # Boundary type for tidal boundary (See boundary types below - Bt sometimes crashes, Bout and Br are prefered)
boundarytype2   = parameters['BoundaryType2'][0]                     # Boundary type for remaining boundaires
discharge       = int(parameters['Discharge_cms'][0])                # River discharge in m3/s
meshes          = parameters['Mesh']                                 # Name of mesh polygons
mesh_scales     = parameters['Scale'].astype(int)                    # Maximum triangle area of meshes
base_scale      = min(mesh_scales) #50000                            # Maximum triangle area in ocean areas
landcovermethod = parameters['LandcoverMethod'][0]                   # Landcover map type

mesh_string = 'Meshes'
for mesh in range(len(meshes)):
    mesh_string = mesh_string +'_' + meshes[mesh] + '_' + str(mesh_scales[mesh]) + 'm2'

##############################################################################
## Polygons to use for mesh refinement should be stored in /Mesh folder
## Mesh names and scales must be listed in the configuration file (Columns 'Mesh' and 'Scale')
## File naming format for polygons should be D_AAAA_mesh.shp where:
## D is the name of your study area
## AAAA is your naming convention for different areas of refinement (rivers, oceans in the example)
## Example is magdalena_rivers_mesh.shp
## For successful mesh generation:
##      1) Polygons should not intersect
##      2) Polygons should contain no holes (we attempt to remove those below)\
##      3) Polygons should not be more complex (fine detail) than the desired resolution (ie you cannot have coarse resolution in very tiny areas)
##      4) Polygons must fall within the model domain
##      5) Polygons will be buffered to ensure transition between areas of different refinement.
##      6) If mesh generation fails, please try simple polygons first, check mesh output, and then add complexity
##

mesh_shapefiles =[]
for mesh in range(len(meshes)):
    mesh_shapefiles = np.concatenate((mesh_shapefiles,[os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(folders[3])
            for f in fnmatch.filter(files,'%s*%s_mesh_%s.shp' %(delta,meshes[mesh],res))]),axis=0)

if len(mesh_shapefiles)==0:
    print('Uniform mesh')
else:
    print('Polygons for mesh refinement are:')
    print(*mesh_shapefiles,sep='\n')

##############################################################################
## Establish folders for output files
model = '%s_%scms_%s_%s%s%s%s' %(delta,discharge,startdate,boundarytype2,boundarytype1,boundarytype1,boundarytype1)
meshname = '%s%s_%s.tsh' %(folders[3],delta,mesh_string)
meshpath = deltapath + 'Scenarios/' + elevationpath.split('/')[-1] + '/' + mesh_string + '/'
modelpath = delta + '/Scenarios/' + elevationpath.split('/')[-1] + '/' + mesh_string + '/' + model + '_' + nowtime + '/'

print('\nModel path will be %s' %(modelpath))
Path('%s' %(meshpath)).mkdir(parents=True, exist_ok=True)
Path('%s' %(local_path + modelpath)).mkdir(parents=True, exist_ok=True)
Path('%s' %(local_path + modelpath + '/plot')).mkdir(parents=True, exist_ok=True)
Path('%s' %(folders[1])).mkdir(parents=True, exist_ok=True)


print('################################# Working directory will be: \n\n%s\n' %(modelpath))

##############################################################################
############################ Make/Get Input Files ############################
##############################################################################

# ######################### Make Manning Roughness File ########################
## Get final landcover map
landcover_file = rasterio.open([os.path.join(dirpath,f)
            for dirpath,dirnames, files in os.walk(elevationpath)
            for f in fnmatch.filter(files,'*landcover_final.tif')][0])
save_profile = landcover_file.profile
landcover = landcover_file.read(1).astype('int')
print('\n##################### Fetching Manning coefficients LUT')
if landcovermethod == 'Copernicus':
    lut = np.genfromtxt(templates_path + 'Copernicus_manningLUT.csv',delimiter = ',')
    #np.savetxt('%s%s_manningLUT.csv'%(folders[0],delta),lut,delimiter=',')
elif landcovermethod == 'WorldCover':
    lut = np.genfromtxt(templates_path + 'WorldCover_manningLUT.csv',delimiter = ',')
    #np.savetxt('%s%s_manningLUT.csv'%(folders[0],delta),lut,delimiter=',')
else:
    lut = np.genfromtxt("%s%s_manningLUT.csv" %(folders[0],delta),delimiter=',')
print('\n##################### Manning roughness coefficients set according to %s_ManningLUT.csv' %(landcovermethod))
manning = lut[landcover]
frict_name = 'Manning'
save_profile['dtype'] = 'float64'
with rasterio.open("%s%s_%s_%s.tif" %(folders[7],delta,frict_name,landcovermethod),'w', **save_profile) as dst:
    dst.write_band(1,manning.astype('float64'))


##############################################################################
################################ Model Extent ################################
##############################################################################

base_triangle  = int((base_scale*2)**0.5)
extentpoly     = gpd.read_file("%s%s_modeldomain.shp" %(folders[6],delta))
ulx,lry,lrx,uly = extentpoly.total_bounds   # Coordiates converted to UTM coordinate system
temp, domain   = getpolygonpoints(extentpoly,ulx,uly,lrx,lry,base_triangle)
perimeter = temp.length
extent = tuple((point.coords[0][0],point.coords[0][1]) for point in (domain.geometry[:]))
extent2 = []
for feat in domain.geometry:
    extent2.append([round(feat.x,5), round(feat.y,5)])
extent2 = np.array(extent2)

try:
    tide_bnd = gpd.read_file("%s%s_tidebnd.shp" %(folders[0],delta))
    print('Tidal boundary is defined by user %s%s_tidebnd.shp' %(folders[0],delta))
except:
    north = np.where((extent2[:,1]==round(uly,5)))
    south = np.where((extent2[:,1]==round(lry,5)))
    west  = np.where((extent2[:,0]==round(ulx,5)))
    east  = np.where((extent2[:,0]==round(lrx,5)))
    other_tag = np.where((extent2[:,1]!=round(uly,5))&(extent2[:,1]!=round(lry,5))&(extent2[:,0]!=round(ulx,5))&(extent2[:,0]!=round(lrx,5)))
else:
    tide_tag = gpd.overlay(extentpoly,tide_bnd,how='intersection')
    other_tag = gpd.overlay(extentpoly,tide_tag,how='difference')



np.savetxt('%s%s_extent_%sm2.csv' %(local_path + modelpath,delta,base_scale), extent,delimiter=',', fmt= '%1.3f') ## USED BY ANUGA MODEL

#bounding_polygon = anuga.read_polygon('../../%s_extent.csv' %(delta))
bounding_polygon = anuga.read_polygon('%s%s_extent_%sm2.csv' %(local_path + modelpath,delta,base_scale))
A = anuga.polygon_area(bounding_polygon) / 1000000.0
print('################################# Bounding Polygon defined by :        %s%s_extent_%sm2.csv' %(modelpath,delta,base_scale))
print('################################# Area of bounding polygon:        %.2f km^2' %(A))
print('################################# Base triangle area for unrefined areas:     %s km^2' %(base_scale/1000000.00))

##############################################################################
##############################################################################
##############################################################################
############################# Build Mesh Polygons ############################
##############################################################################
##############################################################################
##############################################################################

## If polygon complexity is high and vertex spacing is less than length of maximum triangle edge, the mesh generator will fail ##
## There can be NO intersections of polygons
## There should be no holes in the polygons. A  polygon can be completely contained with another, but they must be defined separately
## All polygons must be WITHIN model extent
## The code below will attempt to remove holes and other invalid geometries, but it is not guarenteed
def delete_holes(input_polys):
       crs = input_polys.crs
       exploded = input_polys.explode(index_parts = True)
       interiors = exploded.interiors
       print('There are %s holes' %(len(interiors)))
       holes = gpd.GeoDataFrame()
       holes['geometry'] = None
       holes.crs = crs
       i=0
       for polygon in interiors:
           for interior in range(len(polygon)):
               interior_x = np.array(polygon[interior].coords.xy[0])
               interior_y = np.array(polygon[interior].coords.xy[1])
               newxy = np.column_stack((interior_x,interior_y))
               temppolygon = Polygon(newxy)
               holes.loc[i,'geometry'] = temppolygon
               i=i+1
       holes.geometry = holes.buffer(0.1)
       input_polys = gpd.overlay(input_polys,holes,how='union',keep_geom_type=True)
       #input_polys.geometry = input_polys.buffer(0.1)
       input_polys['dissolve'] = 1
       output_polys = input_polys.dissolve(by='dissolve').reset_index(drop=True)
       print('holes remove')
       return output_polys


interior_regions = []
minlength = (min(mesh_scales)*2)**0.5
maxlength = (max(mesh_scales)*2)**0.5


missing = gpd.GeoDataFrame(columns =['geometry'], crs = EPSG)
for mesh in range(len(mesh_shapefiles)):
    try:
        mesh_scale = mesh_scales[mesh]        #zone_scale = int(parameters['Zone%sScale' %(zone)].astype('float64'))     # Maximum triangle area in Zone A areas
    except:''
    else:
        print('################################# Max triangle area for %s:     %s km^2' %(meshes[mesh],mesh_scale/1000000.00))
        zone_triangle = (mesh_scale*2)**0.5
        zone_polygons = []

        ## ZONE
        zone_polygon = gpd.read_file(mesh_shapefiles[mesh])
        #zone_polygon.geometry = zone_polygon.buffer(2*int(maxlength))
        zone_polygon.geometry = zone_polygon.buffer(-int(minlength))
        zone_no_holes = delete_holes(zone_polygon)
        zone_no_holes = zone_no_holes.explode(index_parts = True).reset_index(drop=True)

        columns = zone_no_holes.columns[:-1]
        zone_no_holes.drop(columns,inplace=True,axis=1)
        zone_filled = zone_no_holes.copy(deep=True)
        points = {'geometry':[Point(zone_no_holes.geometry[0].exterior.coords[0][0], zone_no_holes.geometry[0].exterior.coords[0][1])]}
        points = gpd.GeoDataFrame(points,crs=EPSG)
        rivercsv = []
        for g in range(zone_filled.shape[0]):
            try:
                interior_regions.append([anuga.read_polygon(folders[3]+delta+'_%s%s_points.csv' %(meshes[mesh],g)),mesh_scale])
                points2 = gpd.read_file(folders[3]+delta+'_%s%s_points.csv' %(meshes[mesh],g))
                zone_filled.loc[g].geometry = temp_filled.geometry[0]
                points = pd.concat((points,points_2),axis=0)
                points = points.reset_index(drop=True)
            except:
                tempgdf = gpd.GeoDataFrame(geometry = zone_filled.loc[g])
                #tempgdf.columns = ['geometry']
                temp_filled,temp_points = getpolygonpoints(tempgdf,ulx,uly,lrx,lry,int(zone_triangle))
                if temp_points.shape[0]>3:
                    points_2 = removenearbypoints(temp_points,domain, folders[3],delta,'%s%s' %(meshes[mesh],g),minlength)
                    #rivercsv.append(deltapath+folders[6]+delta+'_river%s_points.csv' %(g))
                    try:
                        interior_regions.append([anuga.read_polygon(folders[3]+delta+'_%s%s_points.csv' %(meshes[mesh],g)),mesh_scale])
                    except: print('BAD POLYGON --> BAD MESH. Excluding this polygon for the mesh generator')
                    else:
                        zone_filled.loc[g].geometry = temp_filled.geometry[0]
                        points = pd.concat((points,points_2),axis=0)
                        points = points.reset_index(drop=True)
            domain = pd.concat((domain,points),axis=0).reset_index(drop=True)
            missing = pd.concat((missing,zone_filled),axis=0)
        points.to_file("%s%s_points_%s.shp" %(folders[3],delta,meshes[mesh]))
## This fills in the buffer zone between refinement areas
## Right now the mesh size of the buffer zone is set at the smallest mesh size
inner_extent = extentpoly.copy(deep=True)
inner_extent.geometry = extentpoly.buffer(-int(minlength)-10)
gaps = gpd.overlay(inner_extent,missing,how='difference')
#gaps = gaps[gaps['Type']=='Polygon']
#gaps.to_file("%s%s%s_gaps.shp" %(deltapath,folders[6],delta))

# print('################################# Max triangle area for gaps:     %s km^2' %(min(mesh_scales)/1000000.00))
# zone_triangle = minlength
# zone_no_holes = delete_holes(gaps)
# zone_no_holes = zone_no_holes.explode(index_parts = True).reset_index(drop=True)
# zone_no_holes.to_file("%s%s_gaps.shp" %(folders[3],delta))

# columns = zone_no_holes.columns[:-1]
# zone_no_holes.drop(columns,inplace=True,axis=1)
# zone_filled = zone_no_holes.copy(deep=True)
# points = {'geometry':[Point(zone_no_holes.geometry[0].exterior.coords[0][0], zone_no_holes.geometry[0].exterior.coords[0][1])]}
# points = gpd.GeoDataFrame(points,crs=EPSG)
# rivercsv = []
# for g in range(zone_filled.shape[0]):
#     tempgdf = gpd.GeoDataFrame(geometry = zone_filled.loc[g])
#     #tempgdf.columns = ['geometry']
#     temp_filled,temp_points = getpolygonpoints(tempgdf,ulx,uly,lrx,lry,int(zone_triangle))
#     if temp_points.shape[0]>200:
#         points_2 = removenearbypoints(temp_points,domain, folders[3],delta,'gaps%s' %(g),minlength)
#         try:
#             interior_regions.append([anuga.read_polygon(folders[3]+delta+'_gaps%s_points.csv' %(g)),min(mesh_scales)])
#         except: print('BAD POLYGON --> BAD MESH. Excluding this polygon for the mesh generator')
#         else:
#             zone_filled.loc[g].geometry = temp_filled.geometry[0]
#             points = points.append(points_2)
#             points = points.reset_index(drop=True)
# points.to_file("%s%s_gaps_points.shp" %(folders[3],delta))
# domain = domain.append(points).reset_index(drop=True)

domain.to_file("%s%s_all_points.shp" %(folders[3],delta))

print('Building Mesh')
geo_reference = Geo_reference(zone=EPSG,
                                datum='wgs84',
                                projection='UTM',
                                false_easting=500000,
                                false_northing=0)
if os.path.isfile("%s%s_tidebnd.shp" %(folders[0],delta)) == True:
        anuga.create_mesh_from_regions(
            bounding_polygon,
            boundary_tags={'Tide': tide_tag,
                           'other': other_tag},
            maximum_triangle_area = base_scale,
            filename =  meshname,
            interior_regions = interior_regions,
            use_cache = False,
            verbose = False,
            poly_geo_reference=geo_reference,
            mesh_geo_reference=geo_reference,
            )
else:
    anuga.create_mesh_from_regions(
        bounding_polygon,
        boundary_tags={'North': north[0],
                       'East': east[0],
                       'South': south[0],
                       'West': west[0]},
        maximum_triangle_area = base_scale,
        filename =  meshname,
        interior_regions = interior_regions,
        use_cache = False,
        verbose = False,
        poly_geo_reference=geo_reference,
        mesh_geo_reference=geo_reference,
        )
mesh_domain = anuga.create_domain_from_file(meshname)

filein = open(templates_path+ 'Run_Model_V0.py')
template = Template(filein.read())
replacements = {'config':'config_%s_%s.csv' %(delta,datetime.now().strftime("%Y%m%d")),
                'SWOT_HYDROLOGY_TOOLBOX':'$SWOT_HYDROLOGY_TOOLBOX',
                'modelpath': modelpath}
print('Mesh created is %s' %(meshname))

runanuga = template.substitute(replacements)
file = open('%s/%s_Run_Model_V0.py' %(local_path + modelpath,delta),'w')
file.write(runanuga)
file.close()

parameters.loc[0,'MeshName'] = mesh_string
parameters.to_csv(local_path + modelpath + '/config_%s_%s.csv' %(delta,datetime.now().strftime("%Y%m%d")))

filein = open(templates_path + 'submit_anuga.sh')
template = Template(filein.read())
replacements = {
                'delta':delta,
                'localpath':hpc_path,
                'modelpath': modelpath,
                'memgb': '50gb',
                'cores': 48,
                'run_script':'%s_Run_Model_V0.py' %(delta)}
print('Job submit file created is submit_job.sh')

runanuga = template.substitute(replacements)
file = open('%s/submit_anuga.sh' %(local_path + modelpath),'w')
file.write(runanuga)
file.close()

fig = plt.figure(figsize=(10,10))
dplotter = animate.Domain_plotter(mesh_domain,plot_dir = local_path + modelpath + '/plot/')
plt.triplot(dplotter.triang,linewidth=0.1)
plt.axis('scaled')
plt.tight_layout()
plt.savefig('%s/plot/mesh%s.png' %(local_path + modelpath,myid))
