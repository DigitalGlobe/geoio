#!/usr/bin/env python

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
