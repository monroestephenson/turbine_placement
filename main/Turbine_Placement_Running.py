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
import geopandas as gpd
import pandas as pd
from qgis.PyQt.QtCore import QSizeF
orange_layer = QgsProject.instance().mapLayersByName('orange_layer')[0]
def border(orange_layer: QgsVectorLayer) -> List[Tuple[float, float]]:
    """
    Creates a list of tuples of the border points of the orange layer.

    Args:
        orange_layer (QgsVectorLayer): The layer whose borders are to be extracted.

    Returns:
        List[Tuple[float, float]]: A list of coordinates representing the border of the layer.
    """
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
    ellipses_series = gpd.GeoSeries(ellipses, crs='EPSG:4326')
    all_ellipses = gpd.GeoDataFrame(geometry=ellipses_series)
    
    # Concatenate the new ellipses GeoDataFrame with the original layer using GeoPandas
    newlayer = gpd.GeoDataFrame(pd.concat([all_ellipses, newlayer], ignore_index=True), crs='EPSG:4326')

    # Save the result
    #newlayer.to_file('path/to/newlayer.shp')
    
    return newlayer
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
            if newlayer.contains(potential_turbine).all():  # Check if any part of newlayer contains the turbine
                if all(not potential_turbine.intersects(packed) for packed in packed_elipses):
                    packed_elipses.append(potential_turbine)
                y += step(semi)  # Ensure y is incremented as intended
            else:
                y += step(semi)  # Increment y even if the turbine is not contained to avoid infinite loop
        x += step(semi)  # Correctly increment x after the inner y loop is completed
    return packed_elipses
def convert_gdf_to_qgsvectorlayer(gdf, layer_name="memory_layer", geom_type="Polygon"):
    # Define the layer with its geometry type and CRS
    crs = gdf.crs.to_wkt()
    layer = QgsVectorLayer(f"{geom_type}?crs={crs}", layer_name, "memory")

    # Get the fields from the GeoDataFrame
    layer_data = layer.dataProvider()

    # Prepare fields, excluding 'geometry' and handling dtype conversions
    fields = QgsFields()
    for name, dtype in zip(gdf.columns, gdf.dtypes):
        if name != 'geometry':
            if dtype.name.startswith('float'):
                qtype = QVariant.Double
            elif dtype.name.startswith('int'):
                qtype = QVariant.Int
            else:
                qtype = QVariant.String
            fields.append(QgsField(name, qtype))

    layer_data.addAttributes(fields.toList())
    layer.updateFields()

    # Create and add features
    for index, row in gdf.iterrows():
        feature = QgsFeature(fields)
        feature.setGeometry(QgsGeometry.fromWkt(row.geometry.to_wkt()))
        # Ensure attributes are added without the 'geometry' column
        attributes = [row[col] for col in gdf.columns if col != 'geometry']
        feature.setAttributes(attributes)
        layer_data.addFeature(feature)

    # Update the layer's extent after adding features
    layer.updateExtents()
    return layer
new_layer=generateNewLayer(orange_layer, (100,100), 23, border(orange_layer))
vl = QgsVectorLayer(new_layer.to_json(),"mygeojson","ogr")
if not vl.isValid():
    print("Layer failed to load!")
else:
    QgsProject.instance().addMapLayer(vl)
    print("Layer added successfully!")