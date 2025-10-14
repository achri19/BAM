#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:13:23 2020

@author: Alchrist
"""


import ee
from ee import batch
from IPython.display import Image
import numpy as np
import geemap
import rasterio


ee.Initialize()
def maskS2clouds(image):
  qa = image.select('QA60');
  ##Bits 10 and 11 are clouds and cirrus, respectively.
  cloudBitMask = 1 << 10
  cirrusBitMask = 1 << 11
  ##Both flags should be set to zero, indicating clear conditions.
  cloudmask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
#  cloudmask = qa.bitwiseAnd(cloudBitMask).eq(0)*(qa.bitwiseAnd(cirrusBitMask).eq(0))
  return image.updateMask(cloudmask).divide(10000)
  #return image.updateMask(qa.bitwiseAnd(cloudBitMask).eq(0)).updateMask(qa.bitwiseAnd(cirrusBitMask).eq(0)).divide(100000)


water = ee.Image('users/alchrist/usgs_ned_watermask_075')
watermask = water.eq(0); 

    
point_geom = list([(-91.34240082169215,29.43580137127807),(-91.34240082169215,29.699476380043233),(-91.51062897110621,29.699476380043233),(-91.51062897110621,29.43580137127807)])
point_geom = ee.Geometry.Polygon(point_geom)
trees = ee.FeatureCollection("users/alchrist/WLD_veg_summer")
eco = ee.FeatureCollection("users/alchrist/EcoGeoMorphicCells_LAND")
labels = ee.Image('users/alchrist/superpixels_labels_land_int2')
#Map = geemap.Map()
#Map.addLayer(watersheds)
months = [' ','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

#Load Sentinel-1 C-band SAR Ground Range collection (log scale, VV, descending)
for i in range(1,13):
    collectionVV = ee.ImageCollection('COPERNICUS/S1_GRD')\
        .filter(ee.Filter.eq('instrumentMode', 'IW'))\
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))\
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))\
        .filterMetadata('resolution_meters', 'equals' , 10)\
        .filterBounds(point_geom)\
        .select('VV');
    #print(collectionVV, 'Collection VV'); 
    
    #Load Sentinel-1 C-band SAR Ground Range collection (log scale, VH, descending)
    collectionVH = ee.ImageCollection('COPERNICUS/S1_GRD')\
        .filter(ee.Filter.eq('instrumentMode', 'IW'))\
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))\
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))\
        .filterMetadata('resolution_meters', 'equals' , 1)\
        .filterBounds(point_geom)\
        .select('VH');
    #print(collectionVH, 'Collection VH');
    
    #Add the SAR images to "layers" in order to display them
    SARVV = collectionVV.filter(ee.Filter.calendarRange(2015,2020,'year'))\
        .filter(ee.Filter.calendarRange(i,i,'month')).mosaic();
    SARVH = collectionVH.filter(ee.Filter.calendarRange(2015,2020,'year'))\
        .filter(ee.Filter.calendarRange(i,i,'month')).mosaic();
    #Apply filter to reduce speckle
    SMOOTHING_RADIUS = 30;
    SARVV_filtered = SARVV.focal_mean(SMOOTHING_RADIUS, 'circle', 'meters');
    SARVH_filtered = SARVH.focal_mean(SMOOTHING_RADIUS, 'circle', 'meters');
    #Map.addLayer(SARVV_filtered, {min:-15,max:0}, 'SAR VV Filtered',0);
    #Map.addLayer(SARVH_filtered, {min:-30,max:-10}, 'SAR VH Filtered',0);

    sar_comp = ee.Image.cat(SARVV_filtered,SARVH_filtered);
    #bands_sar = ['VV','VH'];
    out1 = batch.Export.image.toDrive(sar_comp, description='WLD_S1%s_2014-18_S1-S2_cells' %(months[i]), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000)
    process = batch.Task.start(out1) 
    
    collectionS2 = ee.ImageCollection("COPERNICUS/S2")\
        .filter(ee.Filter.calendarRange(2014,2020,'year'))\
        .filter(ee.Filter.calendarRange(i,i,'month'))\
        .filterBounds(point_geom)\
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 1))\
        .map(maskS2clouds);
    
    comp = collectionS2.mean()
    #ndvi = comp.normalizedDifference(['B8', 'B4']).rename('NDVI');
    #ndvimean = ndvi.addBands(labels).reduceConnectedComponents(ee.Reducer.mean(),'b1')
    # optsar = ee.Image.cat(comp.toFloat(),ndvi.toFloat(),SARVV_filtered.toFloat(),SARVH_filtered.toFloat())
    # optsar = optsar.updateMask(watermask)
    # output = optsar.clip(point_geom)
    out = batch.Export.image.toDrive(comp, description='WLD_S2%s_2014-18_S1-S2_cells' %(months[i]), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000)
    process = batch.Task.start(out) 
    
    withlabels = ee.Image.cat(optsar,labels)
    withlabels = withlabels.updateMask(watermask)
    formeans = optsar.reduceConnectedComponents(reducer = ee.Reducer.mean(), labelBand = 'b1')
    training = optsar.sample(point_geom,scale=10,numPixels=5000)
    clusterer = ee.Clusterer.wekaKMeans(5).train(training)
    result = optsar.cluster(clusterer)
    
    output1 = result.clip(point_geom)
    out1 = batch.Export.image.toDrive(output1, description='WLD_kmeans_%s_2014-18_clouds1_S1-S2_cells' %(months[i]), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000)
    process = batch.Task.start(out1)
  

    
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


