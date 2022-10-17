import json
import cv2
import geopandas as gpd
import rasterio as rio
from osgeo import osr, ogr, gdal
from shapely.geometry import Polygon, mapping, MultiPoint


def spatial_to_pixel(geo_matrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    """
    ul_x= geo_matrix[0]
    ul_y = geo_matrix[3]
    x_dist = geo_matrix[1]
    y_dist = geo_matrix[5]
    pixel = int((x - ul_x) / x_dist)
    line = -int((ul_y - y) / y_dist)
    return pixel, line


def spatial_polygon_to_pixel(raster, geojson, selected_polygon=1):
    
    point_list = MultiPoint(geojson.geometry[selected_polygon][0].exterior.coords)
    
    converted_coords = []
    
    for point in point_list:
        
        x, y = spatial_to_pixel(raster.GetGeoTransform(), point.x, point.y)
        pixel_point = x, y
        converted_coords.append(pixel_point)
    
    return(converted_coords)


def spatial_polygons_to_pixels(raster, geojson):
    
    raster = gdal.Open(raster)
    
    pixel_poly_list = []
    
    for index, polygon in enumerate(geojson.geometry):
        pixel_polygon = spatial_polygon_to_pixel(raster, geojson, selected_polygon=index)
        pixel_poly_list.append(pixel_polygon)
        
    return(pixel_poly_list)


"""
COCO JSON Creation
"""

class CocoJson: 
    def toJSON(self):
        return(json.dumps(self, default=lambda o: o.__dict__, indent = 4))
    
    class CocoImage: 
        pass
    
    class CocoImages: 
        pass
        
    class CocoPolyAnn: 
        pass
    
    class CocoPolyAnns: 
        pass

    
def coco_bbox(polygon):
    """
    Generate a COCO format bounding box from a Polygon.
    """
    
    bounds = polygon.bounds
    top_left_x = bounds[0]
    top_left_y = bounds[1] #lowest y val, cause it's from top down.
    width = bounds[2] - bounds[0]
    height = bounds[3] - bounds[1]
    cc_bbox = [top_left_x, top_left_y, width, height] #https://www.immersivelimit.com/tutorials/create-coco-annotations-from-scratch/#coco-dataset-format
    return(cc_bbox)


def coco_polygon_annotation(pixel_polygon):
    
    annot = CocoJson.CocoPolyAnn()
    annot.segmentation = [item for sublist in pixel_polygon for item in sublist]
    annot.area = Polygon(pixel_polygon).area
    annot.iscrowd = 0
    annot.image_id = 1
    annot.bbox = coco_bbox(Polygon(pixel_polygon))
    # need to read this from the 
    annot.category_id = 1
    annot.id = 1
    
    return(annot)


def coco_polygon_annotations(pixel_poly_list):
    annotations = CocoJson.CocoPolyAnns()
    annotations.annotations = [coco_polygon_annotation(pixel_poly) for pixel_poly in pixel_poly_list]
    return(annotations)


def raster_to_coco(raster_file):
    """
    Generate a COCO format image object.
    """
    
    raster = cv2.imread(raster_file)
    
    image = CocoJson.CocoImage()
    image.license = 1
    image.filename = raster_file
    image.height = raster.shape[0]
    image.width = raster.shape[1]
    image.id = 1

    return(image)


def coco_image_annotations(raster_file_list):
    images = CocoJson.CocoImages()
    images.images = [raster_to_coco(raster_file) for raster_file in raster_file_list]
    return(images)

