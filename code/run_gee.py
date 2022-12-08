#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 15:01:49 2022

@author: alchrist
"""
import sys
delta = sys.argv[1]
extent_file = sys.argv[2]


import ee
from ee import batch
import geopandas as gpd

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

  
resolution  = 10 # Final Output Resolution of all rasters

##############################################################################
########################### Setup Working Directory ##########################
##############################################################################

ee.Initialize()
extent = gpd.read_file(extent_file)
#extent.buffer(5000)
if extent.crs != 'EPSG:4326':
    print('reprojecting to EPSG')
    extent = extent.to_crs('EPSG:4326')

extent.geometry = extent.buffer(0.1,join_style = 2)
exterior_points = extent.geometry[0].exterior.coords[:]
exterior_points = list(map(list, exterior_points))
AOI = ee.Geometry.Polygon(coords=exterior_points)    

print('Getting Sentinel-1 and Sentinel-2 collections')
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

bands = ['NDWI','VV','VH'];

## Process Sentinel-1 Data
VV_ref = collectionVV.max().focal_mean(30, 'circle', 'meters')
VH_ref = collectionVH.max().focal_mean(30, 'circle', 'meters')

# Process Sentinel-2 Data
S2_max = collectionS2.max()    
S2_min = collectionS2.min()
S2_mean = collectionS2.mean()    
ndwi_min = S2_min.select('NDWI')
ndvi_min = S2_min.select('NDVI')   
ndwi_max = S2_min.select('NDWI')
ndvi_max = S2_min.select('NDVI')   
ndwi_mean = S2_min.select('NDWI')
ndvi_mean = S2_min.select('NDVI')   


###############################
###############################
## S2 max, VV, and VH composite
ref_composite = ee.Image.cat(ndwi_max.toFloat(),VV_ref.toFloat(),VH_ref.toFloat()).select(bands).clip(AOI);
## Save composite
out = batch.Export.image.toDrive(ref_composite, folder = 'GEE_drive_composites', description='%s_ndwimax_vvvhmax' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
process = batch.Task.start(out)
print('Saved NDWI max, VV max, and VH max composite to Google Drive at GEE_drive_composites/%s_ndwimax_vvvhmax' %(delta))
## Unsupervised Clustering using S2 max, VV, and VH composite
training = ref_composite.sample(region=AOI,scale=resolution,numPixels=500,tileScale = 16)
classifier = ee.Clusterer.wekaKMeans(2).train(training,bands)
water = ref_composite.select(bands).cluster(classifier);
## Save clusters
out = batch.Export.image.toDrive(water, region=AOI,folder = 'GEE_drive_masks', description='%s_ndwimax_vvvhmax_clustered' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
process = batch.Task.start(out)
print('Saved clusters based on NDWI max, VV max, and VH max to Google Drive at GEE_drive_masks/%s_ndwimax_vvvhmax_clustered' %(delta))

###############################
###############################
## Save S2 mean VV and VH max
ref_composite = ee.Image.cat(ndwi_mean.toFloat(),VV_ref.toFloat(),VH_ref.toFloat()).select(bands).clip(AOI);
## Save composite
out = batch.Export.image.toDrive(ref_composite, folder = 'GEE_drive_composites', description='%s_ndwimean_vvvhmax' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
process = batch.Task.start(out)
print('Saved NDWI mean, VV max, and VH max composite to Google Drive at GEE_drive_composites/%s_ndwimean_vvvhmax' %(delta))
## Unsupervised clustering using s2 mean, vv and vh max
training_clouds = ref_composite.sample(region=AOI,scale=10,numPixels=500,tileScale = 16)
classifier_clouds = ee.Clusterer.wekaKMeans(2).train(training_clouds,bands)
ref_clouds = ref_composite.select(bands).cluster(classifier_clouds);
## Save clusters
out = batch.Export.image.toDrive(ref_clouds, region=AOI, folder = 'GEE_drive_masks', description='%s_ndwimean_vvvhmax_clustered' %(delta), fileFormat = 'GeoTiff',scale = 30,maxPixels=10000000000000 )
process = batch.Task.start(out)
print('Saved clusters based on NDWI mean, VV max, and VH max to Google Drive at GEE_drive_masks/%s_ndwimax_vvvhmean_clustered' %(delta))

###############################
###############################
## Save NDVI and NDWI min
ndwi_ndvi = ee.Image.cat(ndvi_min.toFloat(),ndwi_min.toFloat()).select(['NDWI','NDVI']).clip(AOI);      
ndvi_ndwi_out = batch.Export.image.toDrive(ndwi_ndvi, region=AOI,folder = 'GEE_drive_composites', description='%s_ndwi_ndvi_min' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
process = batch.Task.start(ndvi_ndwi_out)  
print('Saved NDWI min and NDVI min composite to Google Drive at GEE_drive_composites/%s_ndwi_ndvi_min' %(delta))
## Unsupervised Clustering using NDWI and NDVI only       
training = ndwi_ndvi.sample(region=AOI,scale=resolution,numPixels=500,tileScale = 16)
classifier = ee.Clusterer.wekaKMeans(2).train(training,['NDWI','NDVI'])
water2 = ndwi_ndvi.select(['NDWI','NDVI']).cluster(classifier);
## Save clusters
out2 = batch.Export.image.toDrive(water2, region=AOI,folder = 'GEE_drive_masks', description='%s_ndwimin_ndvimin_clustered' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
process = batch.Task.start(out2)
print('Saved clusters based on NDWI min and NDVI min to Google Drive at GEE_drive_masks/%s_ndwimin_ndvimin_clustered' %(delta))
 
  
###############################
###############################
## Save NDVI and NDWI max
ndwi_ndvi = ee.Image.cat(ndvi_max.toFloat(),ndwi_max.toFloat()).select(['NDWI','NDVI']).clip(AOI);      
ndvi_ndwi_out = batch.Export.image.toDrive(ndwi_ndvi, region=AOI,folder = 'GEE_drive_composites', description='%s_ndvi_ndwi_max' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
process = batch.Task.start(ndvi_ndwi_out)  
print('Saved NDWI max and NDVI max composite to Google Drive at GEE_drive_composites/%s_ndvi_ndwi_max' %(delta))
## Unsupervised Clustering using NDWI and NDVI only       
training = ndwi_ndvi.sample(region=AOI,scale=resolution,numPixels=500,tileScale = 16)
classifier = ee.Clusterer.wekaKMeans(2).train(training,['NDWI','NDVI'])
water2 = ndwi_ndvi.select(['NDWI','NDVI']).cluster(classifier);
## Save clusters
out2 = batch.Export.image.toDrive(water2, region=AOI,folder = 'GEE_drive_masks', description='%s_ndwimax_ndvimax_clustered' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
process = batch.Task.start(out2)
print('Saved clusters based on NDWI max and NDVI max to Google Drive at GEE_drive_masks/%s_ndwimax_ndvimax_clustered' %(delta))
  
###############################
###############################
## Save NDVI and NDWI max
ndwi_ndvi = ee.Image.cat(ndvi_mean.toFloat(),ndwi_mean.toFloat()).select(['NDWI','NDVI']).clip(AOI);      
ndvi_ndwi_out = batch.Export.image.toDrive(ndwi_ndvi, region=AOI,folder = 'GEE_drive_composites', description='%s_ndvi_ndwi_mean' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
process = batch.Task.start(ndvi_ndwi_out)  
print('Saved NDWI mean and NDVI mean composite to Google Drive at GEE_drive_composites/%s_ndvi_ndwi_mean' %(delta))
## Unsupervised Clustering using NDWI and NDVI only       
training = ndwi_ndvi.sample(region=AOI,scale=resolution,numPixels=500,tileScale = 16)
classifier = ee.Clusterer.wekaKMeans(2).train(training,['NDWI','NDVI'])
water2 = ndwi_ndvi.select(['NDWI','NDVI']).cluster(classifier);
## Save clusters
out2 = batch.Export.image.toDrive(water2, region=AOI,folder = 'GEE_drive_masks', description='%s_ndwimean_ndvimean_clustered' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
process = batch.Task.start(out2)
print('Saved clusters based on NDWI mean and NDVI mean to Google Drive at GEE_drive_masks/%s_ndwimean_ndvimean_clustered' %(delta))
 



print("We've made several options for water masks, you must download them and pick the best one before continuing")
###############################
# ###############################
# ## Ratios only
# ratios = ee.Image.cat(ndvi_mean.toFloat(),ndwi_mean.toFloat(),VV_ref.toFloat(),VH_ref.toFloat()).select(['NDWI','NDVI','VV','VH']).clip(AOI);      
# ndvi_ndwi_out = batch.Export.image.toDrive(ndwi_ndvi, region=AOI,folder = 'GEE_drive_composites', description='%s_ndvi_ndwi_mean' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
# process = batch.Task.start(ndvi_ndwi_out)  

# ## Unsupervised Clustering using NDWI and NDVI only
# ndwi_ndvi = ee.Image.cat(ndvi._meantoFloat(),ndwi_mean.toFloat()).select(['NDWI','NDVI']).clip(AOI);             
# training = ndwi_ndvi.sample(region=AOI,scale=resolution,numPixels=500,tileScale = 16)
# classifier = ee.Clusterer.wekaKMeans(2).train(training,['NDWI','NDVI'])
# water2 = ndwi_ndvi.select(['NDWI','NDVI']).cluster(classifier);
# ## Save clusters
# out2 = batch.Export.image.toDrive(water2, region=AOI,folder = 'GEE_drive_masks', description='%s_ndwimean_ndvimean_clustered' %(delta), fileFormat = 'GeoTiff',scale = 10,maxPixels=10000000000000 )
# process = batch.Task.start(out2)
          
