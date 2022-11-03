import argparse
import cv2
import geopandas as gpd
import glob
import json
import os.path
import pandas as pd
import rasterio as rio

from pathlib import Path
import os
from itertools import product
from osgeo import osr, ogr, gdal
from rasterio import windows
from shapely.geometry import Polygon, MultiPoint, box


def get_tiles(ds, width=2000, height=2000, map_units=False):

    """
    Defines a set of tiles over a raster layer based on user specified dimensions.

    Args:
        ds: a raster layer that has been read into memory
        width: integer defining the width of tiles
        length: integer defining the length of tiles
        map_units: boolean specifying if width and height at in map units.
    """

    if map_units:
        # Get pixel size
        px, py = ds.transform.a, -ds.transform.e
        width, height = int(width / px + 0.5) , int(height / px + 0.5)

    ncols, nrows = ds.meta['width'], ds.meta['height']

    offsets = product(range(0, ncols, width), range(0, nrows, height))
    big_window = windows.Window(col_off=0, row_off=0, width=ncols, height=nrows)
    for col_off, row_off in  offsets:
        window =windows.Window(col_off=col_off, row_off=row_off, width=width, height=height).intersection(big_window)
        transform = windows.transform(window, ds.transform)
        yield window, transform

def write_raster_tiles(infile, out_path):

    output_filename = "tile_{}-{}.tif"

    with rio.open(infile) as inds:
    
        tile_width, tile_height = 2000, 2000  

        meta = inds.meta.copy()

        for window, transform in get_tiles(inds, tile_width, tile_height, map_units=True):

            meta['transform'] = transform
            meta['width'], meta['height'] = window.width, window.height
            outpath = os.path.join(out_path,output_filename.format(int(window.col_off), int(window.row_off)))
            with rio.open(outpath, 'w', **meta) as outds:
                outds.write(inds.read(window=window))

def get_tile_polygons(raster_tile, geojson, project_crs = "EPSG:3577", filter = True):
    
    """
    Create polygons from a geosjon for an individual raster tile.
    
    Args:
        raster_tile: a file name referring to the raster tile to be loaded
        geojson: a geodataframe with polygons
        
    Returns: 
        tile_polygon: geodataframe with polygons within the raster's extent
    """
    
    # Load raster tile
    raster_tile = rio.open(raster_tile)
    raster_extent = gpd.GeoDataFrame({"id":1,"geometry":[box(*raster_tile.bounds)]}, crs=project_crs)
    geojson = geojson.to_crs(project_crs)
    tile_polygons = geojson.clip(raster_extent)
    # Split multipolygon 
    tile_polygons = tile_polygons.explode(index_parts=False)
    tile_polygons = tile_polygons.reset_index()
    # Filter out zero area polygons
    tile_polygons = tile_polygons[tile_polygons.geometry.area > 0]
    if filter == True:
        tile_polygons = tile_polygons[tile_polygons.geometry.area > 5000]
    tile_polygons = tile_polygons.reset_index()
    
    return(tile_polygons)


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


"""
COCO JSON Creation
"""

class coco_json: 
    def toJSON(self):
        return(json.dumps(self, default=lambda o: o.__dict__, indent = 4))
    
    class coco_image: 
        pass
    
    class coco_images: 
        pass
        
    class coco_poly_ann: 
        pass
    
    class coco_poly_anns: 
        pass
    

def raster_to_coco(raster_file, ind):
    
    """
    Generate a COCO format image object from a raster file.
    """
    
    # TODO: Make this more intelligent.
    try:
        image_extension
    except NameError:
        image_extension = ".png"

    raster = cv2.imread(raster_file)
    # Create a jpg filename and rewrite as a jpg
    raster_name = os.path.splitext(raster_file)[0]
    image_name = f"{raster_name}{image_extension}"
    
    # Write a jpg of the raster tile
    if not os.path.isfile(image_name):
        gdal_command_string = f"gdal_translate -of 'PNG' {raster_file} {image_name}"
        os.system(gdal_command_string)

        # translate_options = gdal.TranslateOptions(format='PNG',
        #                                   outputType=gdal.GDT_Byte,
        #                                   scaleParams=['']
        #                                   )
        # gdal.Translate(destName=image_name, srcDS=raster_file, options=translate_options)


        
    # Create each individual image object
    image = coco_json.coco_image()
    image.license = 1
    image.file_name = os.path.basename(image_name)
    image.height = raster.shape[0]
    image.width = raster.shape[1]
    image.id = ind

    return(image)
            

def coco_image_annotations(raster_file_list):
    
    images = coco_json.coco_images()
    images.images = [raster_to_coco(raster_file, ind) for ind, raster_file in enumerate(raster_file_list)]
    
    return(images)


def spatial_polygon_to_pixel(raster_tile, spatial_polygon):
    
    raster = gdal.Open(raster_tile)
    point_list = MultiPoint(spatial_polygon.exterior.coords)
    converted_coords = []
    
    for point in point_list:
        x, y = spatial_to_pixel(raster.GetGeoTransform(), point.x, point.y)
        pixel_point = x, y
        converted_coords.append(pixel_point)

    return(converted_coords)


def pixel_polygons_for_raster_tiles(raster_file_list, geojson):
    
    tmp_list= []

    for index, file in enumerate(raster_file_list):
        tmp = get_tile_polygons(file, geojson)
        tmp['raster_tile'] = file
        tmp['image_id'] = index
        tmp_list.append(tmp)
        
    pixel_df = pd.concat(tmp_list)
    pixel_df['pixel_polygon'] = pixel_df.apply(lambda row: spatial_polygon_to_pixel(row['raster_tile'], row['geometry']), axis = 1)
    pixel_df['annot_id'] = range(0, 0+len(pixel_df))
    
    return(pixel_df)


def coco_bbox(polygon):

    """
    Generate a COCO format bounding box from a Polygon.
    
    Based on code from:
    #https://www.immersivelimit.com/tutorials/create-coco-annotations-from-scratch/#coco-dataset-format
    """
    
    bounds = polygon.bounds
    top_left_x = bounds[0]
    top_left_y = bounds[1] #lowest y val, cause it's from top down.
    width = bounds[2] - bounds[0]
    height = bounds[3] - bounds[1]
    cc_bbox = [top_left_x, top_left_y, width, height]
    
    return(cc_bbox)


def coco_polygon_annotation(pixel_polygon, image_id, annot_id):
    
    annot = {
        "segmentation":[item for sublist in pixel_polygon for item in sublist],
        "area": Polygon(pixel_polygon).area,
        "iscrowd": 0,
        "image_id": image_id,
        "bbox": coco_bbox(Polygon(pixel_polygon)),
        "category_id": pixel_polygon['class_id'],
        "id": annot_id
            }
    
    return(annot)


def coco_polygon_annotations(polygon_df):
    
    annotations_tmp = []
    for index, row in polygon_df.iterrows():
        annotations_tmp.append(coco_polygon_annotation(row['pixel_polygon'], row['image_id'], row['annot_id']))
        
    return(annotations_tmp)

"""
Create dataset level objects
"""

license_json = {
        "url": "http://creativecommons.org/licenses/by-nc-sa/2.0/",
        "id": 1,
        "name": "Attribution-NonCommercial-ShareAlike License"
    }

info_json = {
        "description": "ePaddocks small test",
        "url": "https://github.com/Sydney-Informatics-Hub/PIPE-3210-Paddock-CV/tree/s4a-dev/data/epaddocks_test",
        "version": "0.1.0",
        "year": 2022,
        "contributor": "Henry Lydecker",
        "date_created": "2022/10/12"
    }

def make_category_object(geojson):
    
    categories = geojson['class'].unique()
    class_ids = geojson['class_id'].unique()
    
    # TODO: Make this less hardcoded
    categories_json = [
        {"supercategory": "other","id": class_ids[0],"name": f"{categories[0]}"},
        {"supercategory": "agriculture","id": class_ids[1],"name": f"{categories[1]}"}
        ]
    
    return(categories_json)


def assemble_coco_json(raster_file_list, geojson, license_json, info_json, categories_json):
    
    pixel_poly_df = pixel_polygons_for_raster_tiles(raster_file_list, geojson)
    
    coco = coco_json()
    coco.images = coco_image_annotations(raster_file_list).images
    coco.annotations = coco_polygon_annotations(pixel_poly_df)
    coco.license = license_json
    coco.categories = categories_json
    coco.info = info_json
    
    return(coco)


#%% Command-line driver

def main(args=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--polygon-file", required=True, default=".", type=Path)
    ap.add_argument("--raster-file", required=True, type=Path)
    ap.add_argument("--tile-size", default = 1000, type=int, help = "Int length in meters of square tiles to generate from raster. Defaults to 1000 meters.")
    ap.add_argument("--tile-dir", required = True, type = Path)
    ap.add_argument("--json-name", default="coco_from_gis.json", type=Path)
    ap.add_argument("--crs", type = str, help = "Specifiy the project crs to use.")
    ap.add_argument("--cleanup", default = False, type = bool, help = "If set to true, will purge *.tif tiles from the directory. Default to false.")
    ap.add_argument("--short-file-name", type = bool, help = "If True, saves a short file name in the COCO for images.")
    args = ap.parse_args(args)

    """
    Read in raster file and save tiles.
    """
    print(f"Creating {args.tile_size} m^2 tiles from {args.raster_file}")
    infile = args.raster_file
    out_path = args.tile_dir
    output_filename = 'tile_{}-{}.tif'

    with rio.open(infile) as inds:
        tile_width, tile_height = args.tile_size, args.tile_size 

        meta = inds.meta.copy()

        for window, transform in get_tiles(inds, tile_width, tile_height, map_units=True):

            meta['transform'] = transform
            meta['width'], meta['height'] = window.width, window.height
            outpath = os.path.join(out_path,output_filename.format(int(window.col_off), int(window.row_off)))
            with rio.open(outpath, 'w', **meta) as outds:
                outds.write(inds.read(window=window))

    # Close the big raster now that we are done with it.
    inds.close()

    # Read raster tiles into a list.
    raster_file_list = []
    for filename in glob.iglob(f'{out_path}/*.tif'):
        raster_file_list.append(filename)

    print(f"{len(raster_file_list)} raster tiles created")
    # Read geojson file.
    geojson = gpd.read_file(args.polygon_file)
    geojson = geojson.to_crs({'init': 'epsg:3577'})

    # Create class_id for category mapping
    geojson['class_id'] = geojson.groupby('class').ngroup()
    categories_json = make_category_object(geojson)

    print("Converting to COCO")
    # We are now ready to make the COCO JSON.
    spatial_coco = assemble_coco_json(raster_file_list, geojson, license_json, info_json, categories_json)

    # Write COCO JSON to file.
    with open (args.json_name, "w") as f:
        f.write(spatial_coco.toJSON())
    print(f"COCO JSON saved to {args.json_name}")

if __name__ == '__main__':
    main()