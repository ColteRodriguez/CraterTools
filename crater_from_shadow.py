# So many packages, some are def not used but being overprepared doesnt hurt i guess
from osgeo import ogr, gdal
from PyQt5.QtWidgets import QMessageBox, QApplication
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

'''
"τ1andτ2aredetermined by selecting samples of bright pixels, 
medium gray pixels, and dark pixels corresponding torocks, 
regolith, and shadows, and computing the gradients between 
the average values of the samples,and they are generally set 
as 60 and 20, respectively" (Li Wu 2018)
'''

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
sampling_values[0], sampling_values[1], sampling_values[2] = 50, 111, 150


# Current predictions are okay, sometimes, the crater rim isnt as dark as a shadow 
# but is just dark enough to not meet the threshold. Consider performing more tests 
# and editing the threshold
#######################################################
# Takes the subsolar azimuth and maximum solar incidence angle as arguments and 
# returns the offset of the brightest reflected pixel
# We define subsolar azimuth as the angle of the sun in top-down (R2)
# counterclockwise of the eastfacing vecor from the point of interest
# /Users/coltenrodriguez/Desktop/fun_cilas_adventure/Color_Sampling.shp
# Return the expected point of greatest reflectance from the illumination condition
def compute_offsets(centroid_x, centroid_y, r, AOI, SA):
    # Two dimensional coordinate offset (CRS)
    xdiff, ydiff = r * np.sin(AOI) * np.cos(SA), r * np.sin(AOI) * np.sin(SA)
    
    # convert offset to pixel coords
    xoff, yoff = int((xdiff - dataset.GetGeoTransform()[0]) / dataset.GetGeoTransform()[1]), int((ydiff - dataset.GetGeoTransform()[3]) / dataset.GetGeoTransform()[5])
    
    return xoff, yoff
#######################################################

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
print(step_light, step_shadow)

# Open the shapefile
shapefile_path = "/Users/coltenrodriguez/Desktop/fun_cilas_adventure/abcdefg/boulders.shp"
shapefile = ogr.Open(shapefile_path, 1)
layer = shapefile.GetLayer()

# Open the raster file
raster_path = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/crater0488/raster/M104827900.tif"
dataset = gdal.Open(raster_path)
band = dataset.GetRasterBand(1)  # Assuming a single band raster

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
            # Nice tripple for loop. loving the O(n^3) runtime
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
