# So many packages, some are def not used but being overprepared doesnt hurt i guess
from qgis.core import (QgsVectorLayer, QgsField, QgsSpatialIndex, QgsFeature,
                       QgsGeometry, QgsPointXY, QgsVectorDataProvider)
from osgeo import ogr, gdal
from qgis.gui import QgsMapTool
from qgis.utils import iface
import matplotlib.pyplot as plt
import qgis.core
from qgis.core import QgsVectorLayer
import rasterio
import geopandas as gpd
from rasterio.features import geometry_mask
import numpy as np
import json
import math


# -> 143 deals with G1/G2 gradient matrices, not important (yet) for crater classification
# establish τ1, τ2, τ3
raster_path = '/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/crater0488/raster/M104827900.tif'
shapefile_path = '/Users/coltenrodriguez/Desktop/fun_cilas_adventure/abcdefg/Color_Sampling.shp'
sampling_values = []

# Open the raster dataset
with rasterio.open(raster_path) as raster_ds:
    # Open the shapefile as a GeoDataFrame
    gdf = gpd.read_file(shapefile_path)
    
    for idx, row in gdf.iterrows():
        # Get the geometry of the current feature
        geometry = row['geometry']
        
        # Mask the raster data using the feature's geometry
        mask = geometry_mask([geometry], out_shape=raster_ds.shape, transform=raster_ds.transform, invert=True)
        
        # Read the pixel values that overlap with the feature
        data = raster_ds.read(1, masked=True)
        values = data[mask]
        
        # Compute and print the average pixel value
        if values.mask.any():
            print(f"Feature {idx}: No data available")
        else:
            avg_value = values.mean()
            sampling_values.append(avg_value)
            print(f"Feature {idx}: Average pixel value = {avg_value:.2f}")
    
gray, light, shadow = sampling_values[0], sampling_values[1], sampling_values[2] + 10.0

#retrieve raster metadata (optional)
# Path to raster json file
path = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/fakeboulders/mapping/raster/M174285569LE.json"
f = open(path)
contents = json.load(f)
#AOI = float(contents["caminfo"]["camstats"]['maximumincidence'])
#SA = 90 - float(contents["caminfo"]["geometry"]["northazimuth"])
AOI = 41.19
SA = 180.0
print(f"MAXIMUM ANGLE OF INCIDENCE, SUBSOLAR AZIMUTH FOR INPUT RASTER: {AOI, SA}")
step_light = [int(math.cos(math.radians(SA))), int(math.sin(math.radians(SA)))]
step_shadow = [-step_light[0], step_light[1]]

# Open the shapefile
shapefile_path = "/Users/coltenrodriguez/Desktop/fun_cilas_adventure/crater_from_shadow2.shp"
shapefile = ogr.Open(shapefile_path, 1)
layer = shapefile.GetLayer()

# Open the raster file
raster_path = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/cpcrater1880/raster/M123784480LE.tif"
dataset = gdal.Open(raster_path)
band = dataset.GetRasterBand(1)  # Assuming a single band raster

# VERY slow
def extract_grad_matrix(band, layer, step_light, step_shadow):
    # Arrays to store suspected craters
    craters = []
    craters_id = []

    # Iterate through each feature in the shapefile
    for feature in layer:
        center_vector = []
        
        # Extract the feature geometry specs, center and min bounding elliplse radius
        geometry = feature.GetGeometryRef()
        centroid = geometry.Centroid()
        
        # The length of the gradient vector
        r = int(math.sqrt(geometry.GetArea() / 3.14) + 5)
        
        # Get the centroid coordinates
        centroid_x, centroid_y = centroid.GetX(), centroid.GetY()
        
        gradient_matrix = []
        gradient_matrix_G2 = []
        for j in range((4 * r)):
            x_cent, y_cent = centroid_x - ((2 * r)), centroid_y + ((2 * r)) - j
            G1_vector = []
            G2_vector = []
            for i in range((4 * r)):
                
                x_light_direction, y_light_direction = x_cent + step_light[0], y_cent + step_light[1]
                x_shadow_direction, y_shadow_direction = x_cent + step_shadow[0], y_cent - step_shadow[1]
                
                pixel_vals = []
                # Get the pixel values
                center_value = band.ReadAsArray(int((x_cent - dataset.GetGeoTransform()[0]) / dataset.GetGeoTransform()[1]), int((y_cent - dataset.GetGeoTransform()[3]) / dataset.GetGeoTransform()[5]), 1, 1)[0, 0]
                pixel_vals.append(center_value)
                light_pixel_value = band.ReadAsArray(int((x_light_direction - dataset.GetGeoTransform()[0]) / dataset.GetGeoTransform()[1]), int((y_light_direction - dataset.GetGeoTransform()[3]) / dataset.GetGeoTransform()[5]), 1, 1)[0, 0]
                pixel_vals.append(light_pixel_value)
                shadow_pixel_value = band.ReadAsArray(int((x_shadow_direction - dataset.GetGeoTransform()[0]) / dataset.GetGeoTransform()[1]), int((y_shadow_direction - dataset.GetGeoTransform()[3]) / dataset.GetGeoTransform()[5]), 1, 1)[0, 0]
                pixel_vals.append(shadow_pixel_value)
                
                # It isnt clear in Li, Wu 2018 whether or not the pixel values are simplified to -1, 0, and 1. I perfrom this 
                # classification here just in case.
                # classify the pixel values
                # Nice quad for loop, loving the O(n^4) runtime
                for i in range(len(pixel_vals)):
                    value = pixel_vals[i]
                    if value <= sampling_values[0]:
                        pixel_vals[i] = -1
                    elif sampling_values[0] < value <= sampling_values[2]:
                        pixel_vals[i] = 0
                    else:
                        pixel_vals[i] = 1
                
                center_vector.append(pixel_vals[1])
                # compute g1, g2, G1, G2
                g1 = pixel_vals[1] - pixel_vals[0]
                g2 = pixel_vals[1] = pixel_vals[2]
                G1 = g1 + g2
                G2 = g1 - g2
                G1_vector.append(G1)
                G2_vector.append(G2)
                x_cent, y_cent = x_cent + step_shadow[0], y_cent + step_shadow[1]
            gradient_matrix.append(G1_vector)
            gradient_matrix_G2.append(G2_vector)
        plt.figure(figsize=(len(gradient_matrix), len(gradient_matrix[0])))
        plt.axis('off')
        plt.imshow(gradient_matrix, cmap = "gray")
        plt.show()
        plt.figure(figsize=(len(gradient_matrix_G2), len(gradient_matrix_G2[0])))
        plt.axis('off')
        plt.imshow(gradient_matrix_G2, cmap = "gray")
        plt.show()
  
######################################################
# Predict misclassified craters
# If neighboring pixels > 15, featuree is a crater 
def walk_from_seed(band, layer):
    for feature in layer:
        
        seed_matrix = []
        feature_mask = []
        
        # Cool geometry stuffs
        geometry = feature.GetGeometryRef()
        centroid = geometry.Centroid()
        x_center, y_center = centroid.GetX(), centroid.GetY()
        r = int(math.sqrt(geometry.GetArea() / 3.14))
        
        # Where to begin walk
        x_start, y_start = x_center - (2 * r), y_center - (2 * r)
        
        # First walk to generate feature mask
        for i in range(8 * r):
            row_vector = []
            feature_row_vector = []
            for j in range(8 * r):
                # read the pixel and add it to the row_vector
                value = band.ReadAsArray(int((x_start - dataset.GetGeoTransform()[0]) / dataset.GetGeoTransform()[1]), int((y_start - dataset.GetGeoTransform()[3]) / dataset.GetGeoTransform()[5]), 1, 1)[0, 0]
                row_vector.append(value)
                
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(x_start, y_start)
                if geometry.Intersects(point):
                    feature_row_vector.append(1)
                else:
                    feature_row_vector.append(0)
                    
                x_start += 1
            # Add the row vector the the seed matrix
            seed_matrix.append(row_vector) 
            feature_mask.append(feature_row_vector)
            y_start += 1
            x_start -= 8 * r
        
        '''
        plt.figure(figsize=(len(seed_matrix), len(seed_matrix[0])))
        plt.axis('off')
        plt.imshow(seed_matrix, cmap = "gray")
        plt.show()
        '''
        
        '''
        It looks to me like Li, Wu 2018 used a depth first search walk from the perimeter pixels 
        to the similar pixels. For now, this alg simply uses the average of all perimeter pixel 
        values as the t3 threshold and compares each pixel in the image to the threshold. It should 
        be noted that this approach causes lots of problems:
        
        1. In instances where features are located close together, the nearby 
            features may be classified as similar pixles and could cause the feature to be mislabeled as a crater
            
        2. +-3 is an arbitrary range applied ot the threshold. For the most accurate classifications,
            This range should change depending on the grayscale variance of the feature and its surroundings, e.g 
            craters are typically located in regions of very low variance so the threshold should be lower -- the 
            opposite is true for boulders -> larger threshold (This should help with the different shades in boulders)
        
        3. DFS avoids double looping (reducing runtime), but is optimized for oop languages, and its worth not re-writing in Java
        '''
        perimeters = []
        similar_count = 0
        for i in range(len(feature_mask)):
            for j in range(len(feature_mask[0])):
                if feature_mask[j][i] == 1 and (feature_mask[j + 1][i] == 0 or feature_mask[j - 1][i] == 0 or feature_mask[j][i + 1] == 0 or feature_mask[j][i - 1] == 0):
                    perimeters.append(seed_matrix[j][i])
        avg = np.mean(perimeters)
        for i in range(len(feature_mask)):
            for j in range(len(feature_mask[0])):
                if avg - 5 < seed_matrix[j][i] < avg + 5 and feature_mask[j][i] != 1:
                    feature_mask[j][i] = 2
                    similar_count += 1
              
        plt.figure(figsize=(len(feature_mask), len(feature_mask[0])))
        plt.axis('off')
        plt.imshow(feature_mask, cmap = "gray")
        plt.show()
        
        image_size = len(feature_mask) * len(feature_mask)
        ratio = (similar_count / image_size) * 100
        if ratio > 10:
            print(f"boulder {feature.GetFID()} is a crater with ratio: {ratio}")
######################################################        
walk_from_seed(band, layer)
