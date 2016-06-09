import imp as _imp
import logging as _logging
import os as _os

# Pull in main modules/classes
from base import GeoImage
from dg import DGImage
import constants
import utils
import downsample

_logger = _logging.getLogger(__name__)

# Pull in optional modules
try:
    _imp.find_module('matplotlib')
    import plotting
except:
    _logger.info('plotting is not available from geoio - matplotlib is missing.')
    pass

# Set version number
_vpath = _os.path.dirname(_os.path.abspath(__file__))
with open(_os.path.join(_vpath,_os.path.pardir,'VERSION.txt')) as _f:
    __version__ = _f.read().strip('\n')