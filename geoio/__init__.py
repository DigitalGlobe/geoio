import imp as _imp
import logging as _logging

# Pull in main modules/classes
from base import GeoImage
from dg import DGImage
import constants
import utils
#import downsample  # numba compile is trigger at import and it necessarily slow

_logger = _logging.getLogger(__name__)

# Pull in optional modules
try:
    _imp.find_module('matplotlib')
    import plotting
except:
    _logger.info('plotting is not available from geoio - matplotlib is missing.')
    pass

# Import version number
from ._version import __version__