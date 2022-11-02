# gis2coco
Convert geospatial datasets of shapefiles and rasters into COCO format datasets of images and JSON.

## User Guide

### GeoJSON to COCO conversion

The script `geojson2coco.py` supports converting a vector shapefile and raster file into a COCO format dataset.

This script assumes that you have one geoJSON with polygon annotations for a geographic extent, and one raster fie that covers this entire geographic extent.

Currently, the supported file formats for vectors are geoJSON and for raster geoTIFF.

To run the conversion, you need to specify where your input data is, as well as where the output should go. 
The script will read in your input data and split the raster into tiles, using `tile-size` to determine how big each tile should be. 

```bash
python scripts/geojson2coco.py \ 
--polygon-file data/example.geojson \
--raster-file data/example.tif \
--tile-size 1000 \
--tile-dir data/example/tiles \
--json-name example_coco.json
```

