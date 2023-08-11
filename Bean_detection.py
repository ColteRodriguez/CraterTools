# Messing around with boulder vs crater geometries. Use crater_from_shadow.py instead

from qgis.core import QgsVectorLayer, QgsProject

# Path to the first shapefile
shapefile1_path = '/Users/coltenrodriguez/Desktop/Stanford_EPSP_23/cpCrater2055_allimgs/CPcrater2055_high_angles/GIS_data/BoulderMap_Mask.shp'
# Path to the second shapefile
shapefile2_path = '//Users/coltenrodriguez/Desktop/Stanford_EPSP_23/cpCrater2055_allimgs/CPcrater2055_high_angles/GIS_data/Bounding_Geometry.shp'

# Load the shapefiles as QgsVectorLayer objects
layer1 = QgsVectorLayer(shapefile1_path, 'Boulder_Map', 'ogr')
layer2 = QgsVectorLayer(shapefile2_path, 'Bounding_geometry', 'ogr')


# Get the feature iterator for both layers
features1 = layer1.getFeatures()
features2 = layer2.getFeatures()


total = 0
# Iterate through each feature in both layers simultaneously
for feature1, feature2 in zip(features1, features2):
    # Get the perimeter attribute value for each feature
    perimeter1 = feature1['Real_Perim']
    perimeter2 = feature2['Ellip_Per']

    # Calculate the difference in perimeters
    ratio = perimeter2/perimeter1

    # Print the difference
    if ratio < 0.95:
        print(feature1['Dist_Rim'], feature2['Dist_Rim'], 'Ratio: {:.2f}'.format(ratio))
        total+=1
        
print("Number of bean shaped boulders: ", total) 
# Clean up and remove the layers from the project
QgsProject.instance().removeMapLayer(layer1)
QgsProject.instance().removeMapLayer(layer2)
