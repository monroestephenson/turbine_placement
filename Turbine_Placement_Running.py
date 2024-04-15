from qgis.core import (
    QgsProject,
    QgsPrintLayout,
    QgsLayoutItemMap,
    QgsPointXY,
    QgsGeometry,
    QgsFeature,
    QgsVectorLayer,
    QgsLayoutExporter
)
from shapely.geometry import Point
from shapely.affinity import scale, rotate
import geopandas as gpd
import pandas as pd
from qgis.PyQt.QtCore import QSizeF
def is_your_point_in_orange_layer(point, orange_layer):
    return orange_layer.intersects(point)
orange_layer = QgsProject.instance().mapLayersByName('orange_layer')[0]
def border(orange_layer): #creates a list of tuples of the border points of the orange layer
    for feature in orange_layer.getFeatures(): 
        geom = feature.geometry()
        if geom.type() == QgsWkbTypes.PolygonGeometry:
            exterior = geom.constGet().boundary()
        border =[]
        for line in exterior:
            vertices = [(vertex.x(), vertex.y()) for vertex in line.vertices()]
            for point in vertices:
                border.append(point)
    return border
def ellipse(center, semi, angle):
    circ = Point(*center).buffer(1)  # Use Point directly since it's already imported
    ell = scale(circ, int(semi[0]), int(semi[1]))
    ellr = rotate(ell, angle)
    elrv = rotate(ell, 90-angle)
    return elrv
def layer_to_gdf(qgs_layer):
    """Convert a QgsVectorLayer to a GeoDataFrame."""
    features = [feature for feature in qgs_layer.getFeatures()]
    gdf = gpd.GeoDataFrame.from_features(features, crs=qgs_layer.crs().toWkt())
    return gdf
def generateNewLayer(orange_layer, semi, angle, border):
    # Ensure the input layer is in the correct CRS
    newlayer = layer_to_gdf(orange_layer)
    # Initialize a list to hold the generated ellipses
    ellipses = []
    
    for point in border:
        elb = ellipse(point, semi, angle)
        # Store the geometry object directly
        ellipses.append(elb)
    
    # Create a GeoSeries from the ellipses with the correct CRS, then convert to GeoDataFrame
    ellipses_series = gpd.GeoSeries(ellipses, crs='EPSG:4326')
    all_ellipses = gpd.GeoDataFrame(geometry=ellipses_series)
    
    # Concatenate the new ellipses GeoDataFrame with the original layer using GeoPandas
    newlayer = gpd.GeoDataFrame(pd.concat([all_ellipses, newlayer], ignore_index=True), crs='EPSG:4326')

    # Save the result
    #newlayer.to_file('path/to/newlayer.shp')
    
    return newlayer
def step(semi):
    return min(semi[0],semi[1])
def place_turbines(newlayer, semi, angle):
    packed_elipses = []
    minx, miny, maxx, maxy = newlayer.total_bounds
    x = minx
    while x < maxx:
        y = miny
        while y < maxy:
            potential_turbine = ellipse((x, y), semi, angle)
            # Adjusted to use .any() for Series boolean context evaluation
            if newlayer.contains(potential_turbine).all():  # Check if any part of newlayer contains the turbine
                if all(not potential_turbine.intersects(packed) for packed in packed_elipses):
                    packed_elipses.append(potential_turbine)
                y += step(semi)  # Ensure y is incremented as intended
            else:
                y += step(semi)  # Increment y even if the turbine is not contained to avoid infinite loop
        x += step(semi)  # Correctly increment x after the inner y loop is completed
    return packed_elipses