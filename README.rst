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

Using the GeoImage object:

.. code:: python

    # Instantiate an image object
    img = geoio.GeoImage('/path/to/imgfile.TIF')

    # Print useful information about the object
    img.files
    img.meta_geoimg

    # Get numpy array
    data = img.get_data()

    # Process data and write to new image
    newdata = data*2
    img.write_img_like_this('/path/to/newfile.TIF',newdata)
    
Using the DGImage object:

.. code:: python

    # Instantiate an image object
    img = geoio.DGImage('/path/to/dgimgfile.TIL')
    # Can also be used directly with a DigitalGlobe TIF file if an XML and/or IMD
    # is available with same name as the TIF file.

    # Print useful information about the object
    img.files
    img.meta_geoimg
    img.meta_dg_quick

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
    data = img.get_data_as_toa_ref()
    
Plotting with the ``geoio.plotting`` functions:

.. code:: python

    # Instantiate an image object
    img = geoio.DGImage('/path/to/dgimgfile.TIL')
    
    # Plot the RGB image
    geoio.plotting.imshow(img.get_data(bands='RGB'))
    
    # Plot the near-infrared false color image
    geoio.plotting.imshow(img.get_data(bands=['N1','G','B']))
    
    # Plotting a histogram of the image bands
    geoio.plotting.hist(img.get_data())
    
    # Plotting a histogram of specific bands
    geoio.plotting.hist(img.get_data(bands='VIS'))
    
