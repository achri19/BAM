#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 10:45:33 2021

@author: Alchrist
Download GEE images from Sentinel-1 and Sentinel -2 for time series analysis
"""
import geopandas as gpd
import ee
import numpy as np
#ee.Authenticate()
def downloadGEE(startyear, endyear,month_interval):
    # Date range and frequency
    years = np.arange(startyear,endyear+1,1)
    months = np.char.zfill(np.arange(1,13,month_interval).astype('str'),2)

    #Initialize GEE and import collections
    ee.Initialize()
    ee.data._use_cloud_api = False
    
    extent = gpd.read_file('/Users/Alchrist/Documents/CurrentProjects/Delta-X/HydroGeomorphic_zones/DeltaX_ROI.shp')
    exterior_points = extent.geometry[0].exterior.coords[:]
    exterior_points = list(map(list, exterior_points))
    AOI = ee.Geometry.Polygon(coords=exterior_points)
    
    training_polygons = ee.FeatureCollection('users/alchrist/clouds')
    watermask = ee.Image('users/alchrist/Sentinel2_DOC_mosaic_sub_EVI_mask_4326')
    
    ## Full collections of Sentinel-1 and Sentinel-2
    collectionVV = ee.ImageCollection('COPERNICUS/S1_GRD')\
        .filter(ee.Filter.eq('instrumentMode', 'IW'))\
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))\
        .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))\
        .filterMetadata('resolution_meters', 'equals' , 10)\
        .filterBounds(AOI)\
        .select('VV');
    collectionVH = ee.ImageCollection('COPERNICUS/S1_GRD')\
        .filter(ee.Filter.eq('instrumentMode', 'IW'))\
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))\
        .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))\
        .filterMetadata('resolution_meters', 'equals' , 10)\
        .filterBounds(AOI)\
        .select('VH');
    lis = ee.List(['15RXN', '15RXP', '15RYN', '15RYP']);
    collectionS2 = ee.ImageCollection("COPERNICUS/S2")\
        .filterBounds(AOI)\
        .filter(ee.Filter.inList('MGRS_TILE',lis))\
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) ##optional
    
    ## Get reference images from April 20-22, which was used to create cloud training polygons
    VV_ref = collectionVV.filterDate('2017-04-20', '2017-04-22').mosaic().focal_mean(50, 'circle', 'meters')
    VH_ref = collectionVH.filterDate("2017-04-20", "2017-04-22").mosaic().focal_mean(50, 'circle', 'meters')
    S2_ref = collectionS2.filterDate('2017-04-20','2017-04-22').mosaic()
    ndvi_ref = S2_ref.normalizedDifference(['B4', 'B8']).rename('NDVI'); # (Red - NIR) / (Red + NIR)
    ndmi_ref = S2_ref.normalizedDifference(['B8','B11']).rename('NDMI'); # (NIR - SWIR) / (NIR + SWIR)
    ref_composite = ee.Image.cat(S2_ref,ndvi_ref,ndmi_ref,VV_ref,VH_ref).clip(AOI);
    bands = ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12','VV','VH','NDVI','NDMI'];
    keepbands = ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12','VV','VH'];
    
    ## Train Random Forest classifier
    training_clouds = ref_composite.select(bands).sampleRegions(training_polygons,['id'],10);
    classifier_clouds = ee.Classifier.smileRandomForest(10).train(training_clouds,'id',bands)
    ref_clouds = ref_composite.select(bands).classify(classifier_clouds);
    
    ## For each yeay, get two-month median mosaics
    #years = ['2017','2018','2019','2020']
    #months = ['01','03','05','07','09','11']
    #months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    
    ## Function to get date of image
    def returndate(image):
        return ee.Date(image).format('YYYY-MM-dd')
    ## Function to get mask cloud classifier
    def getmask(image):
        final_clouds = image.classify(classifier_clouds);
        cloud_mask = final_clouds.eq(5).Or(final_clouds.eq(6)).Or(final_clouds.eq(8)).And(watermask.eq(1))
        return image.updateMask(cloud_mask)
    ## Dummy function to copy image collection and delet placeholder
    def setTemp(image):
        return image.set('Temp',ee.Number(1))
    
    for year in years:
        VV_year = collectionVV.filterDate('%s-01-01' %(year), '%s-12-31' %(year))
        VH_year = collectionVH.filterDate('%s-01-01' %(year), '%s-12-31' %(year))
        OPT_year = collectionS2.filterDate('%s-01-01' %(year), '%s-12-31' %(year))
        for month in months:
            start_date = ee.Date('%s-%s-01' %(year,month))
            end_date = start_date.advance(month_interval,'month')
            print(start_date.format('YYYY-MM-dd').getInfo(),end_date.format('YYYY-MM-dd').getInfo())    
            OPT_month = OPT_year.filterDate(start_date,end_date)
            dates = OPT_month.aggregate_array('system:time_start').map(returndate)
            unique = ee.List(dates).distinct()
            final_collection = ee.ImageCollection.fromImages([OPT_month.mosaic()]);
            final_collection = final_collection.map(setTemp)
            for T in range(unique.length().getInfo()):
                dateOfInterest = unique.get(T)
                def dateDistance(image):
                    return image.set('dateDist',ee.Number(image.get('system:time_start')).subtract(ee.Date(dateOfInterest).millis()).abs());
                #print('date',dateOfInterest.getInfo())
                #start_date = ee.Date(dateOfInterest).advance(-1,'day')
                #end_date = ee.Date(dateOfInterest).advance(1,'day')
                VV = VV_year.map(dateDistance).sort('dateDist',True).limit(4); # Approximately 4 tiles for each date to cover the AOI
                VV_final = VV.mosaic().focal_mean(50, 'circle', 'meters')
                test = VV.first()
                #print('VV',ee.Date(test.get('system:time_start')).format('YYYY-MM-dd').getInfo())
                
                VH = VH_year.map(dateDistance).sort('dateDist',True).limit(4); # Approximately 4 tiles for each date to cover the AOI
                VH_final = VH.mosaic().focal_mean(50, 'circle', 'meters')
                test = VH.first()
                print('VH',ee.Date(test.get('system:time_start')).format('YYYY-MM-dd').getInfo())
                
                Opt = OPT_month.map(dateDistance)
                Opt = Opt.filter(ee.Filter.lt('dateDist',100000000))
                Opt_final = Opt.mosaic()
                test = Opt.first()
                print('Opt',ee.Date(test.get('system:time_start')).format('YYYY-MM-dd').getInfo())
                ndvi_final = Opt_final.normalizedDifference(['B4', 'B8']).rename('NDVI');
                ndmi_final = Opt_final.normalizedDifference(['B8', 'B11']).rename('NDMI');
                
                final_composite = ee.Image.cat(Opt_final,ndvi_final,ndmi_final,VV_final,VH_final);
                final_composite = final_composite.toDouble().clip(AOI)
                temp_collection = ee.ImageCollection.fromImages([final_composite]);
                final_collection = ee.ImageCollection(final_collection.merge(temp_collection))
    
            final_collection = final_collection.filter(ee.Filter.notEquals('Temp',1))
            final = final_collection.map(getmask).select(keepbands)
            final_month = final.reduce(ee.Reducer.median()).clip(AOI)
            task1 = ee.batch.Export.image.toDrive(image=final_month,  # an ee.Image object.
                                                 description='%s%s_10m_%s' %(year,month,month_interval),
                                                 folder='GEE',
                                                 fileNamePrefix='%s%s_10m_%s' %(year,month,month_interval),
                                                 scale=10,
                                                 crs='EPSG:4326',
                                                 fileFormat='GeoTIFF',
                                                 maxPixels = 10000000000000)
            task1.start()


downloadGEE(2015,2020,1)
