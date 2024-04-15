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
import geopandas as gpd

from qgis.PyQt.QtCore import QSizeF

# Create a new project  

orange_layer = OgsProject.instance().layer_by_name('orange')

def create_ellipse(center, width, height):
    ellipse = QgsGeometry.fromWkt(f'POINT({center.x()} {center.y()})').buffer(width/2, 30)
    return ellipse

def is_your_point_in_orange_layer(point, orange_layer):
    return orange_layer.intersects(point)

crs = orange_layer.crs().toWkt()
#writer = QgsVectorLayer('/path/to/new')

for feature in orange_layer.getFeatures():
    geom = feature.geometry()
    if geom.type() == QgsWkbTypes.PolygonGeometry:
        

"""
def border(orange_layer): #creates a list of tuples of the border points of the orange layer
    for feature in orange_layer.getFeatures(): 
        geom = feature.geometry()
        if geom.type() == QgsWkbTypes.PolygonGeometry:
        exterior = geom.constGet().exteriorRing()
        border =[]
        for point in exterior:
            tuple = (point.x(), point.y())
            border.append(tuple)
    return border       
 """

""" 
def ellipse(center, semi, angle):
    circ = shapely.geometry.Point(ellipse[0]).buffer(1)
    ell  = shapely.affinity.scale(circ, int(ellipse[1][0]), int(ellipse[1][1]))
    ellr = shapely.affinity.rotate(ell,ellipse[2])
    elrv = shapely.affinity.rotate(ell,90-ellipse[2])
    return elrv
 """

""" 
def generateNewLayer(orange_layer, semi, angle, border):
    newlayer= orange_layer
    for point in border:
      elb=ellipse(point, semi, angle)
        newlayer = elb+newlayer:
            elb=elb.to_crs('EPSG:4326')
            newlayer = newlayer.to_crs('EPSG:4326')
            newlayer = gpd.pd.concat([elb,newlayer])

    newlayer.to_file('path/to/newlayer.shp')
         """


""" 
def step(semi):
    return min(semi[0],semi[1])
def place_turbines(newlayer, semi, angle, starting):
    packed_elipses = []
    minx, miny, maxx, maxy = newlayer.shape.bounds()
    x= minx
    while x < maxx:
        y= miny
        while y < maxy:
            potential_turbine = ellipse((x,y), semi, angle)
            if newlayer.contains(potential_turbine): #also check that (x,y) \in orange_layer
                if all(not potential_turbine.intersects(packed) for packed in packed_elipses):
                    packed_elipses.append(potential_turbine)
                y+=step(semi)
            x+=step(semi)
    return packed_elipses
 """

""" 
def newshpfile(packed_elipses):
newlayer= empty
for element in packed_elipses:
    element=element.to_crs('EPSG:4326')
    newlayer = newlayer.to_crs('EPSG:4326')
    newlayer = gpd.pd.concat([element,newlayer])
return newlayer
 """


""" 
def main(orange_layer, semi, angle):
    border = border(orange_layer)
    newlayer = generateNewLayer(orange_layer, semi, angle, border)
    packed_elipses = place_turbines(newlayer, semi, angle)
    turbine_Placement = newshpfile(packed_elipses)
    turbine_Placement.to_file('path/to/turbine_Placement.shp')
 """

#main(orange_layer, semi, angle)



######################################################################################

""" 
from shapely.geometry import Point
from shapely.affinity import scale, rotate
import geopandas as gpd
from qgis.core import QgsWkbTypes
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

def get_border(orange_layer):
    border = []
    for feature in orange_layer.getFeatures():
        geom = feature.geometry()
        if geom.type() == QgsWkbTypes.PolygonGeometry:
            exterior = geom.constGet().exteriorRing()
            for point in exterior:
                border.append((point.x(), point.y()))
    return border

def ellipse(center, semi, angle):
    circ = Point(center).buffer(1)
    ell = scale(circ, semi[0], semi[1])
    ellr = rotate(ell, angle)
    elrv = rotate(ellr, 90 - angle)
    return elrv

def generate_new_layer(orange_layer, semi, angle, border):
    newlayer = orange_layer.copy()
    for point in border:
        elb = ellipse(point, semi, angle)
        # Assuming 'elb + newlayer' implies spatial union; this pseudocode is unclear.
        # Here we append 'elb' as a new feature to 'newlayer' instead.
        # Proper GIS operation might require a different approach.
        newlayer = gpd.GeoDataFrame(pd.concat([gpd.GeoSeries(elb), newlayer.geometry], ignore_index=True))
        newlayer.crs = 'EPSG:4326'
    newlayer.to_file('path/to/newlayer.shp')
    return newlayer

def step(semi):
    return min(semi)

def place_turbines(newlayer, semi, angle):
    packed_elipses = []
    bounds = newlayer.total_bounds
    minx, miny, maxx, maxy = bounds
    x = minx
    while x < maxx:
        y = miny
        while y < maxy:
            potential_turbine = ellipse((x, y), semi, angle)
            if newlayer.contains(potential_turbine):  # Spatial check
                if all(not potential_turbine.intersects(packed) for packed in packed_elipses):
                    packed_elipses.append(potential_turbine)
            y += step(semi)
        x += step(semi)
    return packed_elipses

def new_shp_file(packed_elipses):
    newlayer = gpd.GeoDataFrame(crs='EPSG:4326')
    for element in packed_elipses:
        newlayer = gpd.GeoDataFrame(pd.concat([gpd.GeoSeries(element), newlayer.geometry], ignore_index=True))
    newlayer.crs = 'EPSG:4326'
    return newlayer

def main(orange_layer, semi, angle):
    border_points = get_border(orange_layer)
    newlayer = generate_new_layer(orange_layer, semi, angle, border_points)
    packed_elipses = place_turbines(newlayer, semi, angle)
    turbine_placement = new_shp_file(packed_elipses)
    turbine_placement.to_file('path/to/turbine_Placement.shp')




 """