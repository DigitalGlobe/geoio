import unittest

import geoio
import dgsamples

class TestDownsample(unittest.TestCase):
    """Test accuracy of downsampling routines.
    """
    def setUp(self):
        # Setup test gdal object
        self.test_img = dgsamples.wv2_longmont_1k.ms
        self.img = geoio.GeoImage(self.test_img)

    def tearDown(self):
        # Remove gdal image object
        self.img = None

    def test_GeoImage_files_meta_exists(self):
        pass # TBD
        #self.assertIsInstance(self.img.files,tt.bunch.OrderedBunch)


if __name__ == '__main__':
    unittest.main()