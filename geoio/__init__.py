import imp as _imp

# Pull in main modules/classes
from base import GeoImage
from dg import DGImage
import constants

# Pull in optional modules
try:
    _imp.find_module('matplotlib')
    import plotting
except:
    print('plotting is not available from geoio - either skimage or matplotlib '
          'is missing')
    pass
