from typing import List, Tuple
import qgis
import argparse
from qgis.PyQt.QtWidgets import QInputDialog, QLineEdit
import sys
from qgis.core import (
    QgsProject,
    QgsPrintLayout,
    QgsLayoutItemMap,
    QgsPointXY,
    QgsGeometry,
    QgsFeature,
    QgsField, QgsFields, Qgis,
    QgsVectorLayer,
    QgsLayoutExporter
)
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.affinity import scale, rotate
from shapely.geometry import LineString, MultiPoint
import geopandas as gpd
import pandas as pd
from qgis.PyQt.QtCore import QSizeF

#input_layer = QgsProject.instance().mapLayersByName('orange_layer')[0]
def interpolated_points(vertices: QgsVectorLayer, max_distance: float) -> List[Tuple[float, float]]:
    """
    Interpolates points along the lines defined by the vertices of geometries in a QgsVectorLayer,
    ensuring no segment between consecutive points exceeds the specified max_distance.

    Parameters:
        input_layer (QgsVectorLayer): The layer containing line geometries to be processed.
        max_distance (float): The maximum allowed distance between consecutive points.

    Returns:
        List[Tuple[float, float]]: A list of tuples, where each tuple represents the coordinates (x, y)
                                   of interpolated points along the lines of the input_layer.
    """
    interpolated_points = []
    coords = [(vertex.x(), vertex.y()) for vertex in vertices]  # Extract coordinates from QgsPoint
    line = LineString(coords)
    num_segments = int(line.length // max_distance)
    for i in range(num_segments + 1):
        point = line.interpolate(i * max_distance)
        interpolated_points.append((point.x, point.y))
    return interpolated_points
def border(input_layer: QgsVectorLayer, max_distance: float) -> List[Tuple[float, float]]:
    """
    Creates a list of tuples of the border points of the orange layer, including interpolated points
    on straight segments.

    Args:
        input_layer (QgsVectorLayer): The layer whose borders are to be extracted.
        max_distance (float): The maximum distance between points on the border.

    Returns:
        List[Tuple[float, float]]: A list of coordinates representing the detailed border of the layer.
    """
    detailed_border = []
    for feature in input_layer.getFeatures():
        geom = feature.geometry()
        if geom.type() == QgsWkbTypes.PolygonGeometry:
            exterior = geom.constGet().boundary()
            for line in exterior:
                vertices = [vertex for vertex in line.vertices()]  # Extract vertices from the line
                detailed_border.extend([(vertex.x(), vertex.y()) for vertex in vertices])  # Add vertices as tuples
                # Interpolate additional points along straight line segments
                interpolate_points = interpolated_points(vertices, max_distance)
                detailed_border.extend(interpolate_points)
    return detailed_border
def ellipse(center: Tuple[float, float], semi: Tuple[float, float], angle: float) -> Polygon:
    """
    Generates an ellipse from a center point, semi-axis lengths, and rotation angle.

    Args:
        center (Tuple[float, float]): The center of the ellipse.
        semi (Tuple[float, float]): A tuple of semi-major and semi-minor axis lengths.
        angle (float): The rotation angle of the ellipse.

    Returns:
        Polygon: The rotated ellipse as a Shapely polygon.
    """
    circ = Point(*center).buffer(1)  # Use Point directly since it's already imported
    ell = scale(circ, int(semi[0]), int(semi[1]))
    ellr = rotate(ell, angle)
    elrv = rotate(ell, 90-angle)
    return elrv
def layer_to_gdf(qgs_layer: QgsVectorLayer) -> gpd.GeoDataFrame:
    """
    Converts a QgsVectorLayer to a GeoDataFrame.
    """
    if not qgs_layer or not qgs_layer.isValid():
        raise ValueError("Invalid QgsVectorLayer provided.")
    
    features = [feature for feature in qgs_layer.getFeatures()]
    gdf = gpd.GeoDataFrame.from_features(features, crs=qgs_layer.crs().toWkt())
    return gdf

def generateNewLayer(input_layer: QgsVectorLayer, semi: Tuple[float, float], angle: float, border: List[Tuple[float, float]]) -> gpd.GeoDataFrame:
    """
    Generates a new layer consisting of ellipses based on points in the provided border.

    Args:
        input_layer (QgsVectorLayer): The original vector layer.
        semi (Tuple[float, float]): Semi-axis lengths for the ellipses.
        angle (float): Rotation angle for the ellipses.
        border (List[Tuple[float, float]]): Points on the border to place ellipses.

    Returns:
        gpd.GeoDataFrame: The new layer as a GeoDataFrame.
    """
    # Ensure the input layer is in the correct CRS
    newlayer = layer_to_gdf(input_layer)
    # Initialize a list to hold the generated ellipses
    ellipses = []
    
    for point in border:
        elb = ellipse(point, semi, angle)
        # Store the geometry object directly
        ellipses.append(elb)
    
    # Create a GeoSeries from the ellipses with the correct CRS, then convert to GeoDataFrame
    ellipses_series = gpd.GeoSeries(ellipses, crs='EPSG:8857')
    all_ellipses = gpd.GeoDataFrame(geometry=ellipses_series)
    all_ellipses_union = all_ellipses.unary_union
    newlayer_union = newlayer.unary_union
    # Concatenate the new ellipses GeoDataFrame with the original layer using GeoPandas
    #newlayer = gpd.GeoDataFrame(pd.concat([all_ellipses, newlayer], ignore_index=True), crs='EPSG:8857')
    combined_geometries = gpd.GeoSeries([all_ellipses_union, newlayer_union])
    combined=combined_geometries.unary_union
    new_combined_layer = gpd.GeoDataFrame(geometry=combined_geometries, crs='EPSG:8857')

    # Save the result
    #newlayer.to_file('path/to/newlayer.shp')
    
    return new_combined_layer

def place_turbines(newlayer: gpd.GeoDataFrame, semi: Tuple[float, float], angle: float, input_layer: QgsVectorLayer) -> List[Polygon]:
    """
    Places turbines on a new layer as ellipses, ensuring they are packed without overlapping
    and contained within the layer's boundaries.
    """
    packed_elipses = []
    centroids=[]
    minx, miny, maxx, maxy = newlayer.total_bounds
    x = minx
    # Convert the input layer to a GeoDataFrame once before the loop
    input_gdf = layer_to_gdf(input_layer)

    while x < maxx:
        y = miny
        while y < maxy:
            potential_turbine = ellipse((x, y), semi, angle)
            if newlayer.contains(potential_turbine).any() and input_gdf.contains(potential_turbine.centroid).any():
                if all(not potential_turbine.intersects(packed) for packed in packed_elipses):
                    packed_elipses.append(potential_turbine)
                    centroids.append((x,y))
                y += 10
            else:
                y += 10
        x += 10
    return packed_elipses, centroids

def create_layer_from_shapely_polygons(polygons: List[Polygon], layer_name: str = "Polygons", crs: str = "EPSG:8857") -> QgsVectorLayer:
    """
    Creates a new QgsVectorLayer from a list of Shapely Polygon objects.

    Args:
        polygons (list): List of Shapely Polygon objects.
        layer_name (str): Name of the new memory layer.
        crs (str): Coordinate reference system for the new layer.

    Returns:
        QgsVectorLayer: The new memory layer containing the polygons.
    """
    # Create a new memory layer for polygon features
    layer = QgsVectorLayer(f"Polygon?crs={8857}", layer_name, "memory")
    provider = layer.dataProvider()
    
    # Start editing the layer
    layer.startEditing()

    # Create a feature for each Shapely Polygon object
    for shapely_polygon in polygons:
        # Create a new feature
        feature = QgsFeature()
        # Set the geometry of the feature to the Shapely Polygon
        feature.setGeometry(QgsGeometry.fromWkt(shapely_polygon.wkt))
        # Add the feature to the provider
        provider.addFeature(feature)
    
    # Commit changes and update the layer's extent
    layer.commitChanges()
    layer.updateExtents()
    
    return layer

def create_layer_from_points(centroids: List[Tuple[float, float]], layer_name: str, crs: str = "EPSG:8857") -> QgsVectorLayer:
    """
    Creates a new QgsVectorLayer from a list of centroid tuples.

    Args:
        centroids (list): List of tuples representing centroids.
        layer_name (str): Name of the new memory layer.
        crs (str): Coordinate reference system for the new layer.

    Returns:
        QgsVectorLayer: The new memory layer containing the centroids.
    """
    # Convert tuples to Shapely Point objects
    points = [QgsGeometry.fromPointXY(QgsPointXY(*centroid)) for centroid in centroids]
    layer = QgsVectorLayer(f"Point?crs={crs}", layer_name, "memory")
    provider = layer.dataProvider()

    # Start editing the layer
    layer.startEditing()

    # Create a feature for each Point object
    for point in points:
        # Create a new feature
        feature = QgsFeature()
        # Set the geometry of the feature to the Point
        feature.setGeometry(point)
        # Add the feature to the provider
        provider.addFeature(feature)

    # Commit changes and update the layer's extent
    layer.commitChanges()
    layer.updateExtents()

    return layer
def get_user_input() -> Tuple[QgsVectorLayer, Tuple[float, float], float, float]:
    """
    Prompts the user to input parameters required for generating an ellipse on a map layer.

    This function utilizes a series of dialogs to retrieve the following from the user:
    1. The path to the input layer file which is expected to be a named layer within the current QGIS project.
    2. The width of the ellipse in the specified unit (e.g., meters).
    3. The height of the ellipse in the specified unit.
    4. The angle of rotation for the ellipse in degrees.
    5. The maximum allowed distance between interpolated points on the ellipse.

    The function checks if the input layer specified by the user is valid. If the layer is not valid,
    an error message is printed and None is returned.

    Returns:
        tuple: A tuple containing:
              - QgsVectorLayer: The vector layer specified by the user.
              - tuple: A tuple with two floats representing the width and height of the ellipse.
              - float: The angle of the ellipse in degrees.
              - float: The maximum distance allowed between interpolated points.
              If the layer is not valid, returns None.

    Raises:
        RuntimeError: If the user cancels the input dialog, causing the functions to return invalid values.

    Note:
        This function assumes it is used within a QGIS environment where QgsProject and QInputDialog are available.

    Example:
        >>> layer, dimensions, angle, max_dist = get_user_input()
        >>> print(f"Layer: {layer.name()}, Ellipse Dimensions: {dimensions}, Angle: {angle}, Max Distance: {max_dist}")
    """
    # Using a File dialog to get the layer path
    layer_path, _ = QInputDialog.getText(None, "Input Layer", "Enter the path to the input layer file:", QLineEdit.Normal, "")
    #input_layer = QgsVectorLayer(layer_path, "input_layer", "ogr")
    input_layer = QgsProject.instance().mapLayersByName(layer_path)[0]
    if not input_layer.isValid():
        print("Error: The provided layer is not valid.")
        return None

    # Getting ellipse width and height
    #ellipse_width, _ = QInputDialog.getDouble(None, "Ellipse Width", "Enter the width of the ellipse (e.g., 200):", decimals=2)
    #ellipse_height, _ = QInputDialog.getDouble(None, "Ellipse Height", "Enter the height of the ellipse (e.g., 250):", decimals=2)
    rotor_width, _ = QInputDialog.getDouble(None, "Rotor Width", "Enter the rotor width (e.g., 50):", decimals=2)

    # Getting the angle of the ellipse
    ellipse_angle, _ = QInputDialog.getDouble(None, "Ellipse Angle", "Enter the angle of the ellipse in degrees (e.g., 230):", decimals=2)

    # Getting max distance for interpolation
    max_distance = 5#, _ = QInputDialog.getDouble(None, "Max Distance", "Enter the maximum distance allowed between interpolated points (e.g., 5):", decimals=2)

    return input_layer, rotor_width, ellipse_angle, max_distance
def main():
    inputs = get_user_input()
    if inputs is None:
        return
    input_layer, rotor_width, ellipse_angle, max_distance = inputs
    if not input_layer.isValid():
        raise Exception("Layer is no longer valid")
    
    ellipse_width = rotor_width/0.8
    ellipse_height = rotor_width/(2/3)
    ellipse_dimensions = (ellipse_width, ellipse_height)
    
    input_border = border(input_layer, max_distance)
    new_layer = generateNewLayer(input_layer, ellipse_dimensions, 0, input_border)  # Assuming no rotation is needed for the ellipses

    turbines, centroids = place_turbines(new_layer, ellipse_dimensions, 0, input_layer)  # No rotation
    turbines_layer = create_layer_from_shapely_polygons(turbines, "Turbines")

    # Creating layers for the rotor circles and larger rotor circles
    rotor_circles = [Point(centroid).buffer(rotor_width / 2) for centroid in centroids]
    larger_circles = [Point(centroid).buffer(rotor_width * 1.3 / 2) for centroid in centroids]
    rotor_circles_layer = create_layer_from_shapely_polygons(rotor_circles, "Rotor Circles", "EPSG:8857")
    larger_circles_layer = create_layer_from_shapely_polygons(larger_circles, "Extended Rotor Circles", "EPSG:8857")

    centroids_layer = create_layer_from_points(centroids, "Turbine Centroids", "EPSG:8857")


    # Adding layers to the QGIS project
    QgsProject.instance().addMapLayer(turbines_layer)
    QgsProject.instance().addMapLayer(rotor_circles_layer)
    QgsProject.instance().addMapLayer(larger_circles_layer)
    QgsProject.instance().addMapLayer(centroids_layer)


    print("Turbine placement and rotor circles calculated.")

main()
