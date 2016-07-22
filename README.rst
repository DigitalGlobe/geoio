.. image:: https://badge.fury.io/py/geoio.svg
    :target: https://badge.fury.io/py/geoio

Introduction 
============

Geoio provides facilities to easily interact with geospatial
data. The interactions that are supported include data retrieval, spectral
processing, metadata handling, shapefile intersection/extraction, and retrieval
of statistical information. Specific attention has been paid to accessing
DigitalGlobe data and metadata, but the same facilities in this module can be
used to access non-DigitalGlobe data or to build custom processing and
metadata handling for other satellite platforms.

Installation 
============

.. code:: python

    pip install geoio
    
Dependencies will be handled at install if possible.  GDAL is not cleanly
installable via pip so should be handled separately (conda, yum, apt-get, etc.).
The run dependencies are:  ``gdal, xmltodict, pytz, tzwhere, ephem, numpy, tinytools``.  
Additionally, ``dgsamples`` is required for testing and ``matplotlib`` must be
available for the plotting functions to work.

Note for MAC users: if pip fails for ephem, try installing it directly with conda within
the conda virtual environment, i.e.::

   conda install ephem    
   

Imports 
=======

.. code:: python

    import geoio

imports the main classes ``GeoImage`` and ``DGImage`` to the module root.

The ``GeoImage`` class is a relatively thin wrapper around gdal that provides a
pythonic interface for accessing an arbitrary geospatial image format
(generally those supported by gdal plus the DigitalGlobe .TIL format).
Operations supported include reading, writing, chipping, reprojecting, and meta
data access.  The class methods are populated with reasonable defaults and
object interfaces, making image operations less painful so that you can get on
with the important stuff!

The ``DGImage`` class inherits all the capabilites of ``GeoImage`` and adds
DigitalGlobe meta data handling, spectral processing, and band alias data
retrieval.  Therefore, it requires that the input image be a valid DigitalGlobe
image.  This is currently either a .TIL file with the associated meta data files
(.IMD and/or .XML) present in the image directory or a .TIF files with an
identially named .IMD or .XML file.  The metadata is read into an
``OrderedBunch`` object (inherited from the tinytools package) attached to the
instantiated object.

Quick Start
===========

The geoio classes are best used interactively from within ipython where the 
relevant pretty print methods can be triggered.  Meta data information will be 
reutrned regardless of the interpreter, but the readability is currently 
much better in ipython.

The dgsamples repo is used below.  However, all the operations below can be
run on local data by replacing the dgsamples call with a string to the image
location.  From exmaples, instead of typing ```dgsamples.wv2_longmont_1k.ms```,
a local files at ```/path/to/imgfile.TIF``` can be used.

Using the GeoImage object:

.. code:: python

    import dgsamples

    # Instantiate an image object
    img = geoio.GeoImage(dgsamples.bayou_chip.extract_test)  # a TIF file

    # Print useful information about the object
    img.files
    img.meta

    # Get numpy array
    data = img.get_data()

    # Process data and write to new image
    newdata = data*2
    img.write_img_like_this('/path/to/newfile.TIF',newdata)
    
Using the DGImage object:

.. code:: python

    import dgsamples

    # Instantiate an image object
    img = geoio.DGImage(dgsamples.wv2_longmont_1k.ms)    # a TIL file
    # Can also be used directly with a DigitalGlobe TIF file if an XML and/or IMD
    # is available with same name as the TIF file.

    # Print useful information about the object
    img.files
    img.meta

    # Print the full IMD OrderedBunch object
    img.meta_dg.IMD  # tab completeable through the OrderedBunch

    # Return an ImgArr (a numpy array with band meta data handling)
    data = img.get_data()

    # Convert an ImgArr to a pure numpy array
    npdata = np.asarray(data)

    # Return a pure numpy array
    data = img.get_data(meta=False)

    # Get specific bands using aliases - see geoio.constants.DG_BAND_ALIASES for
    # additional aliases.
    data = img.get_data(bands='VIS')

    # Get specific bands using band aliases
    data = img.get_data(bands=['C','Y'])

    # Get image data and convert to TOA reflectance
    data = img.get_data(stype='toa')

Plotting 
========

Plotting with the ``geoio.plotting`` functions:

.. code:: python
    
    import dgsamples

    # Instantiate an image object
    img = geoio.DGImage(dgsamples.wv2_longmont_1k.ms)  # a TIF file
    
    # Plot the RGB image
    geoio.plotting.imshow(img.get_data(bands='RGB'))
    
    # Plot the near-infrared false color image
    geoio.plotting.imshow(img.get_data(bands=['N1','G','B']))
    
    # Plotting a histogram of the image bands
    geoio.plotting.hist(img.get_data())
    
    # Plotting a histogram of specific bands
    geoio.plotting.hist(img.get_data(bands='VIS'))
    
Spatial Resampling
==================

The geoio module has upsampling and downsmapling code that allows the user
to easily resample two images to the same grid for easy multi-image proceesing.

.. code:: python

    import dgsamples
    
    # Import wv3 images
    ms = geoio.DGImage(dgsamples.wv3_longmont_1k.ms)
    swir = geoio.DGImage(dgsamples.wv3_longmont_1k.swir)
    
    # Upsample the swir image
    swir.upsample_like_that(ms,method='nearest')  # default method is bilinear
    
    # Downsample the ms image
    ms.downsample_like_that(swir)  # default method is aggregation
    
    # Or let geoio figure it out
    ms.resample_like_that(swir)
    swir.resample_like_that(ms)

Iterators
=========

The geoio module also provides several iterators to allow easy access to 
yield based portions of a raster file.

.. code:: python

    import dgsamples
    ms = geoio.DGImage(dgsamples.wv2_longmont_1k.ms)
    
    # iterate through vector geometries
    v = dgsamples.wv2_longmont_1k_vectors.poly_geojson_latlon
    [x for x in ms.iter_vector(vector=v,bands='RGB',mask=True)]
    
    # random windows from the image
    [x.shape for x in ms.iter_window_random(win_size=[10,10], no_chips=20)]
    
    # iterate through image with evenly spaced windows based on requested stride
    [x.shape for x in ms.iter_window(win_size=[10,10], stride=[100,100])
    
