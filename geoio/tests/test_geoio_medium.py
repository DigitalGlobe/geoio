import unittest
import os
import shutil
import collections
import inspect
from osgeo import gdalconst
import numpy as np

import geoio.dg
import tinytools as tt
import geoio
import dgsamples

class TestGeoioSpectralFileHandling(unittest.TestCase):
    """Testing for :
    communication between geoio and DGAComp
    """

    def setUp(self):
        # Ensure the previous tests closed down as they were supposed to by
        # removing all spectral files

        # Setup test img

        self.clean_dir = os.path.dirname(dgsamples.wv2_longmont_1k.ms)
        #self.clean_dir = "../data/imagefiles/053792616010_01/053792616010_01_P001_MUL/"
        print("Copying data")
        test_dir = os.path.dirname(os.path.abspath(inspect.stack()[0][1]))
        self.test_dir = test_dir+os.path.sep+"mul_test"
        shutil.copytree(self.clean_dir,self.test_dir)
        self.test_img =  self.test_dir+os.path.sep+\
                         "14JUN20181517-M2AS-053792616010_01_P001.TIL"
        self.img = geoio.dg.DGImage(self.test_img)

    def tearDown(self):
        # Remove all sepctral files
        del self.img
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_create_at_sensor_rad_files(self):
        pre=(self.img.files.rad==None)
        self.img.create_at_sensor_rad_files()
        post=(self.img.files.rad!=None)
        self.assertTrue(pre & post)

    def test_create_toa_ref_files(self):
        pre=(self.img.files.toa==None)
        self.img.create_toa_ref_files()
        post=(self.img.files.toa!=None)
        self.assertTrue(pre & post)

    @unittest.skip("Need to get off IDL before this can be stable.")
    def test_create_dgacomp_ref_files(self):
        pre=(self.img.files.dgacomp==None)
        self.img.create_dgacomp_ref_files()
        post=(self.img.files.dgacomp!=None)
        self.assertTrue(pre & post)

    def test_DGImage_delete_rad_files(self):
        pre=(self.img.files.rad==None)
        self.img.create_at_sensor_rad_files()
        post=(self.img.files.rad!=None)
        self.img.delete_rad_files(test_only=False)
        gone=(self.img.files.rad==None)
        self.assertTrue(pre & post & gone)

    def test_DGImage_delete_toa_ref_files(self):
        pre=(self.img.files.toa==None)
        self.img.create_toa_ref_files()
        post=(self.img.files.toa!=None)
        self.img.delete_toa_ref_files(test_only=False)
        gone=(self.img.files.toa==None)
        self.assertTrue(pre & post & gone)

    @unittest.skip("Need to get off IDL before this can be stable.")
    def test_DGImage_delete_dgacomp_files(self):
        pre=(self.img.files.dgacomp==None)
        self.img.create_dgacomp_ref_files()
        post=(self.img.files.dgacomp!=None)
        self.img.delete_dgacomp_files(test_only=False)
        gone=(self.img.files.dgacomp==None)
        self.assertTrue(pre & post & gone)

    def test_DGImage_delete_all_spectral_files(self):
        pre=((self.img.files.toa==None) &
             (self.img.files.rad==None))
        self.img.create_toa_ref_files()
        self.img.create_at_sensor_rad_files()
        post=((self.img.files.toa!=None) &
              (self.img.files.rad!=None))
        self.img.delete_all_spectral_files(test_only=False)
        gone=((self.img.files.toa==None) &
              (self.img.files.rad==None))
        self.assertTrue(pre & post & gone)

if __name__ == '__main__':
    unittest.main()