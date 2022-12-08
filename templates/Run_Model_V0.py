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
from pathlib import Path
import numpy as np
import scipy
import scipy.interpolate
import pandas as pd
from datetime import datetime
from string import Template
import geopandas as gpd
from shapely.geometry import Polygon, Point
import rasterio
from osgeo import gdal
import shutil
import rasterio
import matplotlib.pyplot as plt

## ANUGA packages
import anuga
from anuga.utilities.plot_utils import Make_Geotif
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.utilities import animate
from anuga import myid, numprocs, finalize, barrier, Inlet_operator
from anuga.parallel.parallel_inlet_operator import Parallel_Inlet_operator
from pathlib import Path


#path = os.path.dirname(os.path.dirname(elevationpath)) +'/'
#print(path)
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000

if os.getcwd().split('/')[1] == 'Volumes':
    path = '/Volumes/FortressL3/ANUGA/GlobalDeltas/'
    build_path = '/Users/Alchrist/Documents/Github/globalcoast/'
    code_path = '/Users/Alchrist/Documents/Github/globalcoast/code/'
    templates_path = build_path + 'templates/'
    config_path = build_path + 'configure_files/'
    cnes_tools_path = '/Volumes/FortressL3/Temp_ANUGA/Tools/swot-hydroloy-toolbox/'
    local_path = path + 'examples/'
    hpc_path = '/scratch_lg/loac_hydro/alchrist/anuga/build/examples/'

else:
    path = '/scratch_lg/loac_hydro/alchrist/anuga/build/'
    build_path = '/scratch_lg/loac_hydro/alchrist/anuga/build/'
    code_path = '/scratch_lg/loac_hydro/alchrist/code/'
    templates_path = build_path + 'templates/'
    config_path = build_path + 'configure_files/'
    cnes_tools_path = '/scratch_lg/loac_hydro/swot-hydrology-toolbox/'
    local_path = path + 'examples/'


sys.path.insert(1, code_path)
from delta_swot_simulations import make_model_output_rasters, get_orbit_files , run_swot_simulator

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
#deltapath = os.path.dirname(os.path.dirname(os.path.dirname(modelpath)))+'/'       # AOI file path
modelpath = '$modelpath'
parameters  = pd.read_csv(local_path + modelpath + '$config',dtype=str,delimiter=',')

delta = parameters['AOI'][0]                                       # Name of the AOI
deltapath = local_path + delta + '/'

##############################################################################
##############################################################################
##############################################################################
## Initiate timer
now = datetime.now()
nowtime = now.strftime('%Y%m%d')
time00 = time.time()
folders = [deltapath+ 'User_Defined_Files/',
           deltapath+ 'tmp/',
           deltapath+ 'Setup_Files/',
           deltapath+ 'Meshes/',
           deltapath+ 'Scenarios/',
           deltapath+ 'Models/',
           deltapath+ 'Setup_Files/Setup_SHP/',
           deltapath+ 'Setup_Files/Setup_RST/',
           deltapath+ 'Setup_Files/Setup_FIG/']

try:
    os.mkdir(folders[1])
except:''
EPSG           = int(parameters['EPSG'][0])                                  # Coordinate system
elevation      = parameters['ElevationName'][0]                     # Name of the elevation file
elevationpath  = folders[4] + elevation
startdate      = (parameters['StartDate'][0])                    # Date to start simulation
enddate        = (parameters['EndDate'][0])                      # Date to start simulation
boundarytype1  = parameters['BoundaryType1'][0]                     # Boundary type for tidal boundary (See boundary types below - Bt sometimes crashes, Bout and Br are prefered)
boundarytype2  = parameters['BoundaryType2'][0]                     # Boundary type for remaining boundaires
B0             = parameters[parameters['North'][0]][0]
B1             = parameters[parameters['East'][0]][0]
B2             = parameters[parameters['South'][0]][0]
B3             = parameters[parameters['West'][0]][0]
discharge      = int(parameters['Discharge_cms'][0])                # River discharge in m3/s
zones          = parameters['Mesh']
zone_scales    = parameters['Scale']
meshes         = parameters['Mesh']
mesh_scales    = parameters['Scale']#.astype(int)
base_scale     = min(mesh_scales)#50000#int(parameters['BaseScale'].astype('float64'))  # Maximum triangle area in ocean areas
tidex          = float(parameters['TideX'][0])#.astype(float)
tidey          = float(parameters['TideY'][0])#.astype(int)
UpstreamX      = float(parameters['UpstreamX'][0])#.astype(int)
UpstreamY      = float(parameters['UpstreamY'][0])#.astype(int)
simu_length    = (datetime.strptime(enddate,'%Y%m%d')-datetime.strptime(startdate,'%Y%m%d')).days
meshname       = parameters['MeshName'][0]
landcovermethod = parameters['LandcoverMethod'][0]

# if myid == 0:
#     shutil.copyfile('config_magdalena_20211203.csv', modelpath + '/' + os.path.basename('config_magdalena_20211203.csv'))

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

if myid == 0:
    print('######################################################################')
    print('##############################################')
    print('################################# Running ANUGA Model')
    print('##############################################')
    print('######################################################################')
    print('')
    print('################################# Area of Interest:                %s' %(delta))
    print('')
    print('################################# Date and Time: %s ' %(now))

    ##############################################################################
    ## Establish folders for output files
    #meshname = '%s%s_%s.tsh' %(folders[3],delta,mesh_string)
    scenario = '%s_%scms_%s_%s_%s%s%s%s' %(delta,discharge,startdate,meshname,boundarytype2,boundarytype1,boundarytype1,boundarytype1)

    print(elevationpath)



    print('')
    print('################################# Working directory will be: \n\n%s\n' %(local_path + modelpath))
    ##############################################################################
    ############################ Make/Get Input Files ############################
    ##############################################################################

    ######################### Make Manning Roughness File ########################
    # Get final landcover map
    # landcover_file = [os.path.join(dirpath,f) for dirpath,dirnames, files in os.walk(folders[7]) for f in fnmatch.filter(files,'*landcover_final.tif')]
    # landcover_file = rasterio.open(landcover_file[0])
    # save_profile = landcover_file.profile
    # save_profile['dtype'] = 'float32'
    # landcover = landcover_file.read(1).astype('int')

    # # # Get Manning roughness lookup file
    # lut = np.genfromtxt("%s%s_manningLUT.csv" %(folders[0],delta),delimiter=',')
    # print('################################# Manning roughness coefficients set according to /%s%s_manningLUT.csv' %(folders[3],delta))
    # manning = lut[landcover]
    # frict_name = 'Manning'
    # with rasterio.open("%s%s_%s.tif" %(folders[7],delta,frict_name),'w', **save_profile) as dst:
    #     dst.write_band(1,manning)


    #bounding_polygon = anuga.read_polygon('%s%s_extent.csv' %(deltapath,delta))
    #A = anuga.polygon_area(bounding_polygon) / 1000000.0
    ##############################################################################
    ##############################################################################
    ##############################################################################
    ## Create domain based on mesh polygons defined above
    geo_reference = Geo_reference(zone=EPSG,
                                    datum='wgs84',
                                    projection='UTM',
                                    false_easting=500000,
                                    false_northing=0)

    domain = anuga.create_domain_from_file('%s%s_%s.tsh' %(folders[3],delta,meshname))

    domain.set_name(scenario + '_%s' %(nowtime)) # Name of sww file
    domain.set_datadir(local_path + modelpath)
    domain.set_low_froude(1)
    domain.set_flow_algorithm('DE1')
    domain.set_minimum_allowed_height(0.02)
    domain.report_water_volume_statistics
    print(domain.statistics())
    domain.store_centroids = True
    ##############################################################################
    ##############################################################################



    ## Bed elevation
    try:
        domain.set_quantity('elevation', filename= ('%s/%s_%s.pts' %(elevationpath,delta,elevation)), location = 'centroids')
    except:
        try:
            anuga.asc2dem('%s/%s_%s.asc' %(folders[1],delta,elevation),use_cache = False, verbose = False)
        except:
            try:
                domain.set_quantity('elevation', filename= ('%s/%s_%s.asc' %(folders[1],delta,elevation)), location = 'centroids')
            except:
                gdal.Translate('%s/%s_%s.asc' %(folders[1],delta,elevation), '%s/%s_%s.tif' %(elevationpath,delta,elevation), format='AAIGrid',noData =-9999)
            #os.system('gdal_translate -of AAIGrid -a_nodata -9999 %s/%s_%s.tif %s/%s_%s.asc ' %(elevationpath,delta,elevation,elevationpath,delta,elevation))
                domain.set_quantity('elevation', filename= ('%s/%s_%s.asc' %(folders[1],delta,elevation)), location = 'centroids')
                print('################################# Elevation make with:     %s_%s.asc' %(delta,elevation))
            else:
                print('################################# Elevation set with:     %s_%s.asc' %(delta,elevation))

        else:
            anuga.dem2pts('%s/%s_%s.dem' %(folders[1],delta,elevation),use_cache = False, verbose = False)
            domain.set_quantity('elevation', filename= ('%s/%s_%s.pts' %(elevationpath,delta,elevation)), location = 'centroids')
            print('################################# Elevation set with:     %s_%s.pts' %(delta,elevation))

    ##############################################################################
    ##############################################################################
    ##############################################################################
    ## Set Initial Conditions

    ## Manning roughness
    try:
        domain.set_quantity('friction', filename= ('%s/%s_manning.asc' %(folders[1],delta)), location='centroids')        # Constant friction
    except:
        #print('translate tif to asc for manning roughness')
        gdal.Translate('%s/%s_manning.asc' %(folders[1],delta),'%s%s_Manning_%s.tif' %(folders[7],delta,landcovermethod),format='AAIGrid',noData=-9999)
        #os.system('gdal_translate -of AAIGrid -a_nodata -9999 %s%s_manning.tif %s/%s_manning.asc ' %(folders[7],delta,elevationpath, delta))
        domain.set_quantity('friction', filename= ('%s/%s_manning.asc' %(folders[1],delta)), location='centroids')        # Constant friction
    print('################################# Friction set with:      %s_manning.asc' %(delta))
    print('')
    print('')
    ## Stage
    domain.set_quantity('stage', -0.1,location='centroids')
    print('################################# Initial stage set with at -0.1m')

    fig = plt.figure(figsize=(10,10))
    dplotter = animate.Domain_plotter(domain,plot_dir = local_path + modelpath + '/plot/')
    plt.triplot(dplotter.triang,linewidth=0.1)
    plt.axis('scaled')
    plt.tight_layout()
    plt.savefig('%s/plot/mesh%s.png' %(local_path + modelpath,myid))

    domain.print_algorithm_parameters()
    sys.stdout.flush()

else:
    domain=None


domain = anuga.distribute(domain,verbose=True)



##############################################################################
##############################################################################
##############################################################################
# Boundary definitions
#############################################################################
## Get Tide Data
## Run pyTMD to get global tidal predictions
## Set downstream boundary conditions
if myid ==0:
    print('')
    print('')
    print('################################# Using tide data from file: %s_tides_lat%s_lon%s_%s.csv' %(delta,np.round(tidey,2),np.round(tidex,2),startdate))
    print('################################# Downstream boundary conditions are set with water level data from pyTMD')
    print('################################# Lat: %s Lon: %s' %(tidey,tidex))

start_tides = (datetime.strptime(startdate,'%Y%m%d') - datetime.strptime('20100101','%Y%m%d')).days * 24 +1
end_tides = start_tides + (simu_length) * 24 +1
boundary_data = np.genfromtxt("%s%s_tides_lat%s_lon%s_%s.csv" %(folders[2],delta,np.round(tidey,2),np.round(tidex,2),20100101),delimiter=',')
stage_time_fun = scipy.interpolate.interp1d(
    boundary_data[:simu_length*24+1, 0],
    boundary_data[start_tides:end_tides, 1])

Bt = anuga.Transmissive_boundary(domain)
Br = anuga.Reflective_boundary(domain)
Bi1 = anuga.Transmissive_stage_zero_momentum_boundary(domain)
Bi2= anuga.Transmissive_momentum_set_stage_boundary(domain,stage_time_fun)
Bi3 = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain,stage_time_fun)
Bout = anuga.Time_boundary(domain,function=lambda t: [stage_time_fun(t),0.0,0.0])
test = pd.DataFrame([[Bt,Br,Bi1,Bi2,Bi3,Bout]],columns=['Bt', 'Br','Bi1', 'Bi2', 'Bi3','Bout'])
boundary1 = test[boundarytype1][0]
boundary2 = test[boundarytype2][0]

domain.set_boundary({'North': test[B0][0], 'East': test[B1][0], 'South': test[B2][0], 'West': test[B3][0]})
if myid == 0:
    fig = plt.figure(figsize = (10,10), facecolor = 'w', edgecolor = 'k')
    plt.plot(boundary_data[:simu_length*24, 0], stage_time_fun(boundary_data[:simu_length*24, 0]), color = 'black')
    plt.xlabel("Simulation Time (seconds)")
    plt.ylabel("Water Level (m)")
    plt.title("Tide Estimations for %s-%s" %(startdate,enddate))
    plt.savefig("%s/%s_tides.png" %(local_path + modelpath + '/plot',delta))

    print('################################# Boundary types are :')
    print('Bt:   Transmissive_boundary')
    print('Bi1:  Transmissive_stage_zero_momentum_boundary')
    print('Bi2:  Transmissive_momentum_set_stage_boundary')
    print('Bi3:  Transmissive_n_momentum_zero_t_momentum_set_stage_boundary')
    print('Bout: Time_boundary')

    print('################################# Boundaries set as: \nNorth: %s\nEast: %s\nSouth: %s\nWest: %s' %(boundary1,boundary2,boundary1,boundary1))

    print('')
## Set upstream discharge boundary
center = (UpstreamX, UpstreamY)
radius = float(min(zone_scales))#5000
fixed_inflow = None
if center[0] !=0 and center[1] !=0:
    #region0 = anuga.Region(domain,center=center,radius=radius,expand_polygon=True)
    line0 = [[center[0],center[1]],[center[0]-5,center[1]]]
    #print(len(region0.indices))
    #if len(region0.indices)>0:
    fixed_inflow = Inlet_operator(domain,line0,Q=discharge,logging=True, description='upstream river', verbose = False)
    if fixed_inflow is not None:
        print('################################# Discharge of %s m^3/s is input at (%s,%s)' %(discharge,center[0],center[1]))
        fixed_inflow.print_statistics()
else: print('there are no rivers in this model')

if myid == 0:
    filein = open(templates_path + '/Merge_V0.py')
    template = Template(filein.read())
    replacements = {'config':'$config','nowtime':nowtime,'modelpath':modelpath}
    runanuga = template.substitute(replacements)
    file = open('%s/%s_Merge_V0.py' %(local_path + modelpath,delta),'w')
    file.write(runanuga)
    file.close()

    filein = open(templates_path + '/PostProcess_V0.py')
    template = Template(filein.read())
    replacements = {'config':'$config',
                    'modelpath':modelpath,
                    'nowtime':nowtime}#,'SWOT_HYDROLOGY_TOOLBOX':'$SWOT_HYDROLOGY_TOOLBOX'}
    runanuga = template.substitute(replacements)
    file = open('%s/%s_PostProcess_V0.py' %(local_path + modelpath,delta),'w')
    file.write(runanuga)
    file.close()

    filein = open(templates_path + '/submit_swot.sh')
    template = Template(filein.read())

    replacements = {'delta':delta +'_swot',
                    'localpath':local_path,
                    'modelpath':local_path + modelpath,
                    'cores':1,
                    'memgb':'250gb',
                    'run_script':'%s_PostProcess_V0.py' %(delta)}
    print('Job submit file created is submit_job.sh')
    runanuga = template.substitute(replacements)
    file = open('%s/submit_swot.sh' %(local_path + modelpath),'w')
    file.write(runanuga)
    file.close()
##############################################################################
##############################################################################
##############################################################################
## Start simulation
if myid==0:
    startnow = datetime.now()
    print('################################# Starting Simulation at %s' %(startnow.strftime("%d/%m/%Y %H:%M:%S")))
    print('################################# Simulation will run from %s to %s' %(startdate, enddate))
    print('################################# Simulation length will be %s days, %s seconds' %(simu_length,simu_length*24*60*60))
    print('################################# SWW file will be: ' + scenario + '_%s' %(nowtime))# Evolve system through time
else: ''#print('Domain distributed', + myid)

sys.stdout.flush()

#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=3600, duration=(simu_length*24*60*60)):
    if myid==0:
        domain.print_timestepping_statistics()
        sys.stdout.flush()

endnow = datetime.now()

print('################################# Completed Simulation at %s' %(endnow.strftime("%d/%m/%Y %H:%M:%S")))


domain.sww_merge(verbose=True)

anuga.finalize()

if myid ==0:
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
