import imp as _imp
import logging

# Pull in main modules/classes
from base import GeoImage
from dg import DGImage
import constants

logger = logging.getLogger(__name__)

# Pull in optional modules
try:
    _imp.find_module('matplotlib')
    import plotting
except:
    logger.info('plotting is not available from geoio - matplotlib is missing.')
    pass
