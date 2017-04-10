Crop a GeoTIFF
--------------

This example shows how to crop out a subset of a `DGImage`.

This is a slightly modified version of the sample code provided in [issue 11](https://github.com/DigitalGlobe/geoio/issues/11).

The `window` parameter to `get_data` specifies the subset of the image to pull out.  The format is `[LEFT, TOP, WIDTH, HEIGHT]`.  The `return_location` parameter is necessary in order to update the `geo_transform`.

```Python
import geoio
import dgsamples
import pdb

# Create a DGImage
img = geoio.DGImage(dgsamples.wv2_longmont_1k.ms)

cropped = img.get_data(window=[10, 10, 20, 20], return_location=True)

dtype = geoio.constants.DICT_NP_TO_GDAL[cropped[0].dtype]

geoio.base.create_geo_image('./cropped.tif',
                            cropped[0],
                            'GTiff',
                            cropped[1]['geo_transform'],
                            img.meta.projection_string,
                            dtype)

```