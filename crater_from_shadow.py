from osgeo import ogr, gdal
from PyQt5.QtWidgets import QMessageBox
import qgis.core
import numpy as np
import json

# Comments
# Current predictions are okay, sometimes, the crater rim isnt as dark as a shadow 
# but is just dark enough to not meet the threshold. Consider performing more tests 
# and editing the threshold
#######################################################
# Takes the subsolar azimuth and maximum solar incidence angle as arguments and 
# returns the offset of the brightest reflected pixel
# We define subsolar azimuth as the angle of the sun in top-down (R2)
# counterclockwise of the eastfacing vecor from the point of interest
def compute_offsets(dataset, r, AOI, SA):
    # Two dimensional coordinate offset (CRS)
    xdiff, ydiff = r * np.sin(AOI) * np.cos(SA), r * np.sin(AOI) * np.sin(SA)
    
    # convert offset to pixel coords
    xoff, yoff = int((xdiff - dataset.GetGeoTransform()[0]) / dataset.GetGeoTransform()[1]), int((ydiff - dataset.GetGeoTransform()[3]) / dataset.GetGeoTransform()[5])
    
    return xoff, yoff
#######################################################


#retrieve raster metadata
path = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/cpCrater2055_allimgs/CPcrater2055_high_angles/GIS_data/M161855819LE.json"
f = open(path)
contents = json.load(f)
AOI = float(contents["caminfo"]["camstats"]['maximumincidence'])
SA = 90 - float(contents["caminfo"]["geometry"]["northazimuth"])
print(f"MAXIMUM ANGLE OF INCIDENCE, SUBSOLAR AZIMUTH FOR INPUT RASTER: {AOI, SA}")


# Open the shapefile
shapefile_path = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/cpCrater2055_allimgs/CPcrater2055_high_angles/GIS_data/shadow_tests.shp"
shapefile = ogr.Open(shapefile_path, 1)
layer = shapefile.GetLayer()

# Open the raster file
raster_path = "/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/cpCrater2055_allimgs/CPcrater2055_high_angles/GIS_data/M161855819LE_final.tif"
dataset = gdal.Open(raster_path)
band = dataset.GetRasterBand(1)  # Assuming a single band raster

# Arrays to store suspected craters
craters = []
craters_id = []

# Iterate through each feature in the shapefile
for feature in layer:
    
    # Extract the feature geometry specs
    geometry = feature.GetGeometryRef()
    centroid = geometry.Centroid()
    r = feature["r"]

    # Get the pixel coordinates of the centroid
    x, y = centroid.GetX(), centroid.GetY()
    pixel_x = int((x - dataset.GetGeoTransform()[0]) / dataset.GetGeoTransform()[1])
    pixel_y = int((y - dataset.GetGeoTransform()[3]) / dataset.GetGeoTransform()[5])
    
    # Compute the appropriate offstet for illumination condition 
    xoff, yoff = compute_offsets(dataset, r, AOI, SA)
    
    # Read the pixel value
    centroid_pixel_value = band.ReadAsArray(pixel_x, pixel_y, 1, 1)[0, 0]

    # shadow_pixel_value
    # Get the pixel coordinates of the shadowpt
    x, y = centroid.GetX() - r, centroid.GetY()
    pixel_x = int((x - dataset.GetGeoTransform()[0]) / dataset.GetGeoTransform()[1])
    pixel_y = int((y - dataset.GetGeoTransform()[3]) / dataset.GetGeoTransform()[5])
    
    # Read the pixel value
    shadow_pixel_value = band.ReadAsArray(pixel_x, pixel_y, 1, 1)[0, 0]

    # computer the brightness ratio
    ratio = centroid_pixel_value/shadow_pixel_value
    
    # For debugging
    print(f"ID: {feature['id']}. Centroid Pixel Brightness {centroid_pixel_value}, shadow pixel brightness {shadow_pixel_value}, xoff, yoff {xoff, yoff}, brightness ratio {ratio}")
    
    # Threshold allows for some fluctuations in brightness ratio 
    if ratio <= 1.10:
        print("Feature ", feature["id"], "is a boulder with brightness ratio ", centroid_pixel_value/shadow_pixel_value)
        craters.append(feature)
        craters_id.append(feature["id"])
        print("")
        
        
# UI to delete detected craters
msg_box = QMessageBox()
msg_box.setIcon(QMessageBox.Question)
msg_box.setWindowTitle("Confirmation")
msg_box.setText(f"The following features will be deleted: {craters_id}. Continue?")
msg_box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
msg_box.setDefaultButton(QMessageBox.No)

# Display the notification dialog
user_choice = msg_box.exec_()

# Execute deletions
if user_choice == QMessageBox.Yes:
    for feature in craters:
        layer.DeleteFeature(feature.GetFID())

# Cleanup
shapefile = None
dataset = None
