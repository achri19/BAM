#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 09:08:05 2021

@author: Alchrist
"""

import ee
from ee import batch
import geemap
Map = geemap.Map()
from IPython.display import Image
import numpy as np
import geemap
import rasterio
from pathlib import Path
import pandas as pd
import os
import geopandas as gpd
import fnmatch
from osgeo import gdal, osr
import math

path = '/Users/Alchrist/Documents/CurrentProjects/ANUGA/GlobalDeltas/'
configfile = path + 'configure_files/Deltas_SetupParameters_Master_V2.csv'
parameter_file = pd.read_csv(configfile)
deltas = pd.read_csv(configfile)['AOI']

def maskS2clouds(image):
    qa = image.select('QA60');
    ##Bits 10 and 11 are clouds and cirrus, respectively.
    cloudBitMask = 1 << 10
    cirrusBitMask = 1 << 11
    ##Both flags should be set to zero, indicating clear conditions.
    cloudmask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
    #cloudmask = qa.bitwiseAnd(cloudBitMask).eq(0)*(qa.bitwiseAnd(cirrusBitMask).eq(0))
    return image.updateMask(cloudmask)#.divide(10000)
    #return image.updateMask(qa.bitwiseAnd(cloudBitMask).eq(0)).updateMask(qa.bitwiseAnd(cirrusBitMask).eq(0)).divide(100000)  
def addNDWI(image):
  ndwi = image.normalizedDifference(['B8', 'B3']).rename('NDWI').toFloat();
  return image.addBands(ndwi)    
def addNDVI(image):
  #ndvi = image.normalizedDifference(['B8', 'B3']).rename('NDVI').toFloat();
  ndvi = (image.select(['B8']).subtract(image.select(['B3']))).divide(image.select(['B8']).add(image.select(['B3']))).rename('NDVI')
  return image.addBands(ndvi)
def addMNDWI(image):
  mndwi = image.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  return image.addBands(mndwi)
#deltas= ['burdekin','magdalena','gabon','guayas','fly']

for delta in deltas[:]:
    delta = delta.lower()
    print(delta)
    deltapath = path + 'examples/%s/' %(delta)
    resolution = 10 # Final Output Resolution of all rasters
    parameters = parameter_file[parameter_file['AOI']==delta.capitalize()].reset_index()

    ##############################################################################
    ########################### Setup Working Directory ##########################
    ##############################################################################

    while os.path.isdir(deltapath) ==False:
        deltapath = input('The directory %s does not exist, what is the working directory:  ' %(deltapath))
    else:
        print('\n##################### The working directory set as %s ' %(deltapath))
    folders = [deltapath + 'User_Defined_Files/',
            deltapath + 'tmp/',
            deltapath + 'Setup_Files/',
            deltapath + 'Meshes/',
            deltapath + 'Scenarios/',
            deltapath + 'Models/',
            deltapath + 'Setup_Files/Setup_SHP/',
            deltapath + 'Setup_Files/Setup_RST/',
            deltapath + 'Setup_Files/Setup_FIG/']

    ee.Initialize()
    gee = [os.path.join(dirpath,f)
        for dirpath,dirnames, files in os.walk(folders[0])
        for f in fnmatch.filter(files,'%s*input.shp' %(delta))]
    extent = gpd.read_file(gee[0])
    #extent.buffer(5000)
    extent.crs = 'EPSG:4326'
    extent.geometry = extent.buffer(0.1,join_style = 2)
#    extent = extent.to_crs('EPSG:4326')
    exterior_points = extent.geometry[0].exterior.coords[:]
    exterior_points = list(map(list, exterior_points))
    #print(exterior_points)
    AOI = ee.Geometry.Polygon(coords=exterior_points)    
    #AOI.projection().getInfo()
    #AOI = AOI2.buffer(0)
    
    # rivers = gpd.read_file(deltapath + folders[1] + delta + '_watermask.shp')    
    # rivers = rivers.to_crs('EPSG:4326').explode()
    # notrivers = gpd.overlay(extent,rivers,how='difference').explode()
    
    # exterior_points = rivers.geometry[0][0].exterior.coords[:]
    # exterior_points = list(map(list, exterior_points))
    # rivers = ee.Geometry.Polygon(coords = exterior_points)
    
    # exterior_points = notrivers.geometry[0][0].exterior.coords[:]
    # exterior_points = list(map(list, exterior_points))
    # notrivers = ee.Geometry.Polygon(coords = exterior_points)
    collectionS2 = ee.ImageCollection("COPERNICUS/S2")\
        .filter(ee.Filter.calendarRange(2014,2021,'year'))\
        .filter(ee.Filter.calendarRange(1,12,'month'))\
        .filterBounds(AOI)\
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 5))\
        .map(addNDWI)\
        .map(addMNDWI)\
        .map(addNDVI)\
        .map(maskS2clouds)  ;  
    collectionVV = ee.ImageCollection('COPERNICUS/S1_GRD')\
        .filter(ee.Filter.eq('instrumentMode', 'IW'))\
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))\
        .filterMetadata('resolution_meters', 'equals' , 10)\
        .filterBounds(AOI)\
        .select('VV');
    collectionVH = ee.ImageCollection('COPERNICUS/S1_GRD')\
        .filter(ee.Filter.eq('instrumentMode', 'IW'))\
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))\
        .filterMetadata('resolution_meters', 'equals' , 10)\
        .filterBounds(AOI)\
        .select('VH');

    VV_ref = collectionVV.max().focal_mean(30, 'circle', 'meters')
    VH_ref = collectionVH.max().focal_mean(30, 'circle', 'meters')
    bands = ['NDWI','VV','VH'];

    # USE NDWI MAX (over classifies water)
    S2_max = collectionS2.max()    
    ref_composite = ee.Image.cat(S2_max.toFloat(),VV_ref.toFloat(),VH_ref.toFloat()).select(bands).clip(AOI);
    out = batch.Export.image.toDrive(ref_composite, folder = 'GEE_drive/composites', description='%s_ndwimax_vvvhmax' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
    process = batch.Task.start(out)
    
    #Unsupervised Clustering
    training = ref_composite.sample(region=AOI,scale=resolution,numPixels=500,tileScale = 16)
    classifier = ee.Clusterer.wekaKMeans(2).train(training,bands)
    water = ref_composite.select(bands).cluster(classifier);

    out = batch.Export.image.toDrive(water, region=AOI,folder = 'GEE_drive/masks', description='%s_ndwimax_vvvhmax_clustered' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
    process = batch.Task.start(out)
 

    S2_min = collectionS2.min()
    
    ndwi = S2_min.select('NDWI')
    # ndwi_water = ndwi.updateMask(ndwi.gt(0.1)).clip(AOI)
    ndwi_out = batch.Export.image.toDrive(ndwi, region=AOI,folder = 'GEE_drive/composites', description='%s_ndwimin' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
    process = batch.Task.start(ndwi_out)
    
    ndvi = S2_min.select('NDVI')   
    # ndvi_water = ndvi.updateMask(ndvi.gt(0.1)).clip(AOI)
    ndvi_out = batch.Export.image.toDrive(ndvi, region=AOI,folder = 'GEE_drive/composites', description='%s_ndvimin' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
    process = batch.Task.start(ndvi_out)  
    
    
    ndwi_ndvi = ee.Image.cat(ndvi.toFloat(),ndwi.toFloat()).select(['NDWI','NDVI']).clip(AOI);             
    training = ndwi_ndvi.sample(region=AOI,scale=resolution,numPixels=500,tileScale = 16)
    classifier = ee.Clusterer.wekaKMeans(2).train(training,['NDWI','NDVI'])
    water2 = ndwi_ndvi.select(['NDWI','NDVI']).cluster(classifier);
    
    out2 = batch.Export.image.toDrive(water2, region=AOI,folder = 'GEE_drive/masks', description='%s_ndwimin_ndvimin_clustered' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
    process = batch.Task.start(out2)
    # ## Supervised Random Forest classifier
    # points1 = ee.FeatureCollection.randomPoints(rivers,500, 0, 10);
    # points0 = ee.FeatureCollection.randomPoints(notrivers,2, 0, 10);
    # def addClass1(feature):
    #     return feature.set('class1', 1);
    # def addClass0(feature):
    #     return feature.set('class1', 0);
    # points1 = points1.map(addClass1)
    # points0 = points0.map(addClass0)
    # points = points1.merge(points0)
    # training_clouds = ref_composite.sampleRegions(collection=points,properties=['class1'],scale=10)
    # classifier_clouds = ee.Classifier.smileRandomForest(1).train(training_clouds,'class1',bands)
    # ref_clouds = ref_composite.select(bands).classify(classifier_clouds);
    # out = batch.Export.image.toDrive(ref_clouds, folder = 'GEE_drive', description='%s_ndwimax_vvvhmax_classified' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
    # process = batch.Task.start(out)
    
    


    ## USE NDWI MIN (under classifies water)
    S2_ref = collectionS2.mean()    
    ref_composite = ee.Image.cat(S2_ref.toFloat(),VV_ref.toFloat(),VH_ref.toFloat()).select(bands).clip(AOI);
    # out = batch.Export.image.toDrive(ref_composite, folder = 'GEE_drive/composites', description='%s_ndwimean_vvvhmax' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
    # process = batch.Task.start(out)
    
    
    ## Supervised Random Forest classifier
    # points1 = ee.FeatureCollection.randomPoints(rivers,500, 0, 10);
    # points0 = ee.FeatureCollection.randomPoints(notrivers,2, 0, 10);
    # def addClass1(feature):
    #     return feature.set('class1', 1);
    # def addClass0(feature):
    #     return feature.set('class1', 0);
    # points1 = points1.map(addClass1)
    # points0 = points0.map(addClass0)
    # points = points1.merge(points0)
    # training_clouds = ref_composite.sampleRegions(collection=points,properties=['class1'],scale=10)
    # classifier_clouds = ee.Classifier.smileRandomForest(1).train(training_clouds,'class1',bands)
    # ref_clouds = ref_composite.select(bands).classify(classifier_clouds);
    # out = batch.Export.image.toDrive(ref_clouds, folder = 'GEE_drive', description='%s_ndwimean_vvvhmax_classified' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
    # process = batch.Task.start(out)
    
    
    #Unsupervised Clustering
    training_clouds = ref_composite.sample(region=AOI,scale=10,numPixels=500,tileScale = 16)
    classifier_clouds = ee.Clusterer.wekaKMeans(2).train(training_clouds,bands)
    ref_clouds = ref_composite.select(bands).cluster(classifier_clouds);
     
    out = batch.Export.image.toDrive(ref_clouds, region=AOI, folder = 'GEE_drive/masks', description='%s_ndwimean_vvvhmax_clustered' %(delta), fileFormat = 'GeoTiff',scale = 30,maxPixels=10000000000000 )
    process = batch.Task.start(out)
    




    # def filterSpeckles_VV(image):
    #     vv = image.select('VV')
    #     vv_smoothed = vv.focal_median(30,'circle','meters').rename('VV_Filtered') 
    #     return image.addBands(vv_smoothed) 
    # def filterSpeckles_VH(image):
    #     vh = image.select('VH')
    #     vh_smoothed = vh.focal_median(30,'circle','meters').rename('VH_Filtered') 
    #     return image.addBands(vh_smoothed) 

    
    # def lin_to_db(image):
    #     ##Convert backscatter from linear to dB.
    #     bandNames = image.bandNames().remove('angle')
    #     db = ee.Image.constant(10).multiply(image.select(bandNames).log10()).rename(bandNames)
    #     return image.addBands(db, None, True)
    # def db_to_lin(image):
    #     ##Convert backscatter from dB to linear.
    #     bandNames = image.bandNames().remove('angle')
    #     lin = ee.Image.constant(10).pow(image.select(bandNames).divide(10)).rename(bandNames)
    #     return image.addBands(lin, None, True)
    # def maskAngLT452(image):
    #     ##mask out angles >= 45.23993
    #     ang = image.select(['angle'])
    #     return image.updateMask(ang.lt(45.23993)).set('system:time_start', image.get('system:time_start'))
    # def maskAngGT30(image):
    #     ##mask out angles <= 30.63993
    #     ang = image.select(['angle'])
    #     return image.updateMask(ang.gt(35)).set('system:time_start', image.get('system:time_start'))    
    #     #return image.updateMask(ang.gt(30.63993)).set('system:time_start', image.get('system:time_start'))    
    # def maskEdge(image):
    #     ##Remove edges
    #     mask = image.select(0).unitScale(-25, 5).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100)
    #     return image.updateMask(mask.select(0)).set('system:time_start', image.get('system:time_start'))     
    # def f_mask_edges(image):
    #     db_img = lin_to_db(image)
    #     output = maskAngGT30(db_img)
    #     output = maskAngLT452(output)
    #     #output = image.clip(output)
    #     output = maskEdge(output)
    #     output = db_to_lin(output)
    #     return output.set('system:time_start', image.get('system:time_start'))
    # def f_mask_edges2(image):
    #     img = image.select(['VV'])
    #     return image.updateMask(img.gt(-20))
        
    
    # SARVV_filtered = collectionVV.focal_mean(30, 'circle', 'meters');
    # SARVH_filtered = collectionVH.focal_mean(30, 'circle', 'meters');
    #inputs = ee.Image.cat(SARVV_filtered.toFloat(),SARVH_filtered.toFloat()).clip(AOI)

    #.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))\










    with rasterio.open('/Users/Alchrist/Downloads/congo_ndwi(3).tif') as src:
        ndwi = src.read(1)
        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                shapes(ndwi, mask=None, transform=src.transform)))
    geoms = list(results)
    import geopandas as gp
    gpd_polygonized_raster  = gp.GeoDataFrame.from_features(geoms)

    dst_layername = "POLYGONIZED_STUFF"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource( dst_layername + ".shp" )
    dst_layer = dst_ds.CreateLayer(dst_layername, srs = None )
    
    gdal.Polygonize( ndwi, None, dst_layer, -1, [], callback=None )
    
    #MeanImageClip = MeanImage.clip(point_geom)

    #Unsupervised Kmean Clustering
    
    #optsar = ee.Image.cat(S2,SARVV_filtered,SARVH_filtered)
    #sar = SARVV_filtered.addBand(SARVH_filtered)
    bands_opt_sar = ['VV','VH','NDVI'];
    training = optsar.select(bands_opt_sar).sample(point_geom,scale=10,numPixels=5000)
    clusterer = ee.Clusterer.wekaKMeans(5).train(training)
    result = optsar.cluster(clusterer)
    
    output = result.clip(point_geom)
    out = batch.Export.image.toDrive(sar_comp.select(bands_sar), description='WLD_kmeans_%s_2014-18_S1-S2_cells' %(months[i]), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000)
    process = batch.Task.start(out)
    
    Map.addLayer(result.randomVisualizer(), {}, 'classfied')
    
    
#Classify with Sentinel-2 Optical data
opt_comp = ee.Image.cat(S2)
bands_opt= ['B1','B2','B3','B4','B5','B6','B7']
training_opt = opt_comp.select(bands_opt).sampleRegions(trees,['id'],scale=10)
#Train the classifier
classifier_opt = ee.Classifier.randomForest().train(training_opt,'id',bands_opt)
#Run the Classification
classified_opt = composite.select(bands_opt).classify(classifier_opt)
#Display the Classification
#Create a confusion matrix representing resubstitution accuracy.
#print('RF-S2 error matrix: ', classifier_opt.confusionMatrix());
#print('RF-S2 accuracy: ', classifier_opt.confusionMatrix().accuracy());

#Classify with both S1 SAR and S2 Optical data
opt_sar = ee.Image.cat(composite, SARVV_filtered,SARVH_filtered);
bands_opt_sar = ['VV','VH','NDVI','B1','B2','B3','B4','B5','B6','B7'];
training_opt_sar = opt_sar.select(bands_opt_sar).sampleRegions(trees,['id'],scale=10)
#Train the classifier
classifier_opt_sar = ee.Classifier.randomForest().train(training_opt_sar, 'id',bands_opt_sar)
#Run the Classification
classifiedboth = opt_sar.select(bands_opt_sar).classify(classifier_opt_sar);
#Create a confusion matrix representing resubstitution accuracy.
#print('RF-opt-sar error matrix: ', classifier_opt_sar.confusionMatrix());
#print('RF-opt-sar accuracy: ', classifier_opt_sar.confusionMatrix().accuracy());


# output1 = classifiedboth.clip(point_geom)
# out = batch.Export.image.toDrive(output1, description='WLD_classified_summer2014-18_S1-S2', fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000)
# process = batch.Task.start(out)


## Unsupervised Classification


