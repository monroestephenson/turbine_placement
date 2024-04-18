from typing import List, Tuple
import qgis
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
from shapely.geometry import LineString
import geopandas as gpd
import pandas as pd
from qgis.PyQt.QtCore import QSizeF
orange_layer = QgsProject.instance().mapLayersByName('orange_layer')[0]
def interpolate_points(vertices, max_distance):
    """
    Interpolate points along a line defined by QgsPoint vertices if the distance between points is greater than max_distance.
    """
    interpolated_points = []
    coords = [(vertex.x(), vertex.y()) for vertex in vertices]  # Extract coordinates from QgsPoint
    line = LineString(coords)
    num_segments = int(line.length // max_distance)
    for i in range(num_segments + 1):
        point = line.interpolate(i * max_distance)
        interpolated_points.append((point.x, point.y))
    return interpolated_points
def border(orange_layer: QgsVectorLayer, max_distance: float) -> List[Tuple[float, float]]:
    """
    Creates a list of tuples of the border points of the orange layer, including interpolated points
    on straight segments.

    Args:
        orange_layer (QgsVectorLayer): The layer whose borders are to be extracted.
        max_distance (float): The maximum distance between points on the border.

    Returns:
        List[Tuple[float, float]]: A list of coordinates representing the detailed border of the layer.
    """
    detailed_border = []
    for feature in orange_layer.getFeatures():
        geom = feature.geometry()
        if geom.type() == QgsWkbTypes.PolygonGeometry:
            exterior = geom.constGet().boundary()
            for line in exterior:
                vertices = [vertex for vertex in line.vertices()]  # Extract vertices from the line
                detailed_border.extend([(vertex.x(), vertex.y()) for vertex in vertices])  # Add vertices as tuples
                # Interpolate additional points along straight line segments
                interpolated_points = interpolate_points(vertices, max_distance)
                detailed_border.extend(interpolated_points)
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

    Args:
        qgs_layer (QgsVectorLayer): The QGIS vector layer to convert.

    Returns:
        gpd.GeoDataFrame: The converted GeoDataFrame.
    """
    features = [feature for feature in qgs_layer.getFeatures()]
    gdf = gpd.GeoDataFrame.from_features(features, crs=qgs_layer.crs().toWkt())
    return gdf
def generateNewLayer(orange_layer: QgsVectorLayer, semi: Tuple[float, float], angle: float, border: List[Tuple[float, float]]) -> gpd.GeoDataFrame:
    """
    Generates a new layer consisting of ellipses based on points in the provided border.

    Args:
        orange_layer (QgsVectorLayer): The original vector layer.
        semi (Tuple[float, float]): Semi-axis lengths for the ellipses.
        angle (float): Rotation angle for the ellipses.
        border (List[Tuple[float, float]]): Points on the border to place ellipses.

    Returns:
        gpd.GeoDataFrame: The new layer as a GeoDataFrame.
    """
    # Ensure the input layer is in the correct CRS
    newlayer = layer_to_gdf(orange_layer)
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
def step(semi):
    """
    Calculates the step size for iterating over area based on the smaller semi-axis.

    Args:
        semi (Tuple[float, float]): Semi-axis lengths.

    Returns:
        float: The smaller value of the semi-axis lengths.
    """
    return min(semi[0],semi[1])
def place_turbines(newlayer: gpd.GeoDataFrame, semi: Tuple[float, float], angle: float) -> List[Polygon]:
    """
    Attempts to place turbines on a new layer as ellipses, ensuring they are packed without overlapping
    and entirely contained within the layer's boundaries.

    Args:
        newlayer (gpd.GeoDataFrame): The GeoDataFrame on which turbines are to be placed.
        semi (Tuple[float, float]): Semi-axis lengths of the ellipse representing the turbine.
        angle (float): Rotation angle for the ellipses.

    Returns:
        List[Polygon]: A list of successfully placed turbine ellipses.
    """
    packed_elipses = []
    minx, miny, maxx, maxy = newlayer.total_bounds
    x = minx
    while x < maxx:
        y = miny
        while y < maxy:
            potential_turbine = ellipse((x, y), semi, angle)
            # Adjusted to use .any() for Series boolean context evaluation
            if newlayer.contains(potential_turbine).any():  # Check if any part of newlayer contains the turbine
                if all(not potential_turbine.intersects(packed) for packed in packed_elipses):
                    packed_elipses.append(potential_turbine)
                y += 2  # Ensure y is incremented as intended
            else:
                y += 2  # Increment y even if the turbine is not contained to avoid infinite loop
        x += step(semi)  # Correctly increment x after the inner y loop is completed
    return packed_elipses
new_layer=generateNewLayer(orange_layer, (100,150), 23, border(orange_layer,5))

def create_layer_from_shapely_polygons(polygons, layer_name="Polygons", crs="EPSG:4326"):
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
    layer = QgsVectorLayer(f"Polygon?crs={crs}", layer_name, "memory")
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

# Example usage:
# Replace your_list_of_shapely_polygons with the actual list of polygons you have
# Your list of polygons
#for polygon in new_layer:
#    print(polygon.wkt)
# Add the layer to the QGIS interface
QgsProject.instance().addMapLayer(qgs_layer)

vl = QgsVectorLayer(new_layer.to_json(),"mygeojson","ogr")
if not vl.isValid():
    print("Layer failed to load!")
else:
    QgsProject.instance().addMapLayer(vl)
    print("Layer added successfully!")
turbine_placement=place_turbines(generateNewLayer(orange_layer, (100,150), 23, border(orange_layer,5)),(100,150),23)
qgs_layer = create_layer_from_shapely_polygons(turbine_placement)
QgsProject.instance().addMapLayer(qgs_layer)
print(turbine_placement)