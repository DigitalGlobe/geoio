import unittest
import os
import collections
from osgeo import gdalconst
import numpy as np
import warnings

import geoio.dg
import tinytools as tt
import geoio
from geoio import constants as const
import dgsamples

class TestGeoioEnv(unittest.TestCase):
    """Testing for :
    geoio environment
    """

    def test_for_gdal(self):
        #tt.cmd_line.exec_cmd(['gdalinfo','--version'])
        # This is basically testing the the function runs without errors
        #self.assertTrue(True)
        self.assertTrue(tt.cmd_line.is_tool('gdalinfo'))

    def test_for_const(self):
        self.assertEqual(const.DG_ABSCAL_GAIN['WV02_P'], [0.960])
        self.assertEqual(const.DG_ABSCAL_GAIN['WV03_MULTI'][2], 0.907)

    def test_const_dict_np_to_gdal(self):
        a = np.array([1.23],dtype='float32')
        self.assertEqual(const.DICT_NP_TO_GDAL[a.dtype], gdalconst.GDT_Float32)

    def test_const_dict_gdal_to_np(self):
        a = np.array([1.23],dtype='float32')
        self.assertEqual(const.DICT_GDAL_TO_NP[gdalconst.GDT_Float32], a.dtype)

class TestGeoImage(unittest.TestCase):
    """Testing for :
    GeoImage class
    """
    # Run defined method setUp to setup the test environment
    def setUp(self):
        self.test_img = dgsamples.wv2_longmont_1k.ms
        self.img = geoio.GeoImage(self.test_img)

    def tearDown(self):
        # Remove gdal image object
        self.img = None
        # Find temp files from test directory and remove
        tmp_files=tt.files.search(os.path.dirname(os.path.realpath(__file__)),
                        ["*temp*","*tmp*"],case_sensitive=False)

        print(tmp_files)
        [os.remove(f) for f in tmp_files]

    def test_GeoImage_files_meta_exists(self):
        self.assertIsInstance(self.img.files,tt.bunch.OrderedBunch)

    def test_GeoImage_meta_exists(self):
        self.assertIsInstance(self.img.meta_geoimg,tt.bunch.OrderedBunch)

    def test_GeoImage_shape(self):
        self.assertEqual(self.img.shape,(8,500,501))

    def test_GeoImage_resolution(self):
        self.assertEqual(self.img.resolution[0],2)
        self.assertEqual(self.img.resolution[1],2)

    def test_GeoImage_get_data(self):
        a = self.img.get_data()
        self.assertTrue(np.array_equal(a[:,100,100],[472,379,563,690,364,696,447,734]))
        self.assertTrue(np.array_equal(a[:,200,200],[575,439,659,784,402,742,449,706]))
        self.assertTrue(np.array_equal(a[:,300,300],[313,192,253,279,161,376,308,571]))
        self.assertTrue(np.array_equal(a[:,400,400],[322,199,258,283,143,356,287,495]))

    def test_GeoImage_get_data_dtypes(self):
        a = self.img.get_data()
        self.assertEqual(const.DICT_NP_TO_GDAL[a.dtype],
                         self.img.meta_geoimg.data_type)

    def test_GeoImage_write_img_like_this(self):
        a = (self.img.get_data()*0.01).astype('float32')
        back = self.img.write_img_like_this("tmp.tif",a,return_obj=True)
        b = back.get_data()
        self.assertIsNone(np.testing.assert_array_almost_equal(a,b,decimal=6))

class TestGeoImage_getdata(unittest.TestCase):
    # Run defined method setUp to setup the test environment
    def setUp(self):
        self.test_img = dgsamples.wv2_longmont_1k.ms
        self.img = geoio.GeoImage(self.test_img)

    def tearDown(self):
        # Remove gdal image object
        self.img = None
        # Find temp files from test directory and remove
        tmp_files = tt.files.search(os.path.dirname(os.path.realpath(__file__)),
                                    ["*temp*", "*tmp*"], case_sensitive=False)

        print(tmp_files)
        [os.remove(f) for f in tmp_files]

    def test_get_data_no_options(self):
        d = self.img.get_data()
        self.assertEqual(d.shape,(8,501,500))

    def test_get_data_bands_num(self):
        d = self.img.get_data(bands=[7,3,1])
        self.assertEqual(d.shape,(3,501,500))

    def test_get_data_window_centered(self):
        d = self.img.get_data(window=[10,10,21,23])
        self.assertEqual(d.shape,(8,23,21))

    def test_get_data_window_edge(self):
        d = self.img.get_data(window=[499,500,1,1])
        known = np.array(
            [[[354]],[[226]],[[304]],[[341]],[[176]],[[430]],[[312]],[[518]]])
        self.assertTrue(np.array_equal(d,known))

class TestGeoImage_iter_vector(unittest.TestCase):
    # Run defined method setUp to setup the test environment
    def setUp(self):
        self.test_img = dgsamples.wv2_longmont_1k.ms
        self.img = geoio.GeoImage(self.test_img)
        self.vec = dgsamples.wv2_longmont_1k_vectors.poly_geojson_latlon
        self.badvec = dgsamples.bayou_vectors.poly

    def tearDown(self):
        # Remove gdal image object
        self.img = None
        # Find temp files from test directory and remove
        tmp_files = tt.files.search(os.path.dirname(os.path.realpath(__file__)),
                                    ["*temp*", "*tmp*"], case_sensitive=False)

        print(tmp_files)
        [os.remove(f) for f in tmp_files]

    def test_get_data_noOptions(self):
        for x in self.img.iter_vector(vector=self.vec):
            self.assertIsInstance(x,np.ndarray)

    def test_get_data_noOptions_badVec(self):
        for x in self.img.iter_vector(vector=self.badvec):
            self.assertIsNone(x)

    def test_get_data_propTrue(self):
        for x in self.img.iter_vector(vector=self.vec,properties=True):
            self.assertIsInstance(x,tuple)
            self.assertIsInstance(x[0], np.ndarray)
            self.assertIsInstance(x[1], dict)

    def test_get_data_propTrue_badVec(self):
        for x in self.img.iter_vector(vector=self.badvec,properties=True):
            self.assertIsInstance(x,tuple)
            self.assertIsNone(x[0])
            self.assertIsInstance(x[1],dict)

    def test_get_data_propStr(self):
        for x in self.img.iter_vector(vector=self.vec,properties='testfloat'):
            self.assertIsInstance(x,tuple)
            self.assertIsInstance(x[0], np.ndarray)
            self.assertIsInstance(x[1], dict)
            self.assertIsInstance(x[1]['testfloat'], float)

    def test_get_data_propList(self):
        for x in self.img.iter_vector(vector=self.vec,
                                      properties=['testfloat', 'teststr']):
            self.assertIsInstance(x[0], np.ndarray)
            self.assertIsInstance(x[1], dict)
            self.assertIsInstance(x[1]['teststr'], str)
            self.assertIsInstance(x[1]['testfloat'], float)
            self.assertEqual(len(x[1]), 2)

    def test_get_data_propListBad(self):
        with warnings.catch_warnings(record=True) as w:
            for x in self.img.iter_vector(vector=self.vec,
                                          properties=['testfloatBad']):
                self.assertIsInstance(x[0], np.ndarray)
                self.assertIsNone(x[1])
                self.assertTrue(len(w) == 1)

    def test_get_data_propListGoodBad(self):
        for x in self.img.iter_vector(vector=self.vec,
                                      properties=['testfloatBad','teststr']):
            self.assertIsInstance(x, tuple)
            self.assertIsInstance(x[0], np.ndarray)
            self.assertIsInstance(x[1], dict)
            self.assertIsInstance(x[1]['teststr'], str)
            self.assertRaises('warn.warning')
            self.assertEqual(len(x[1]),1)

    def test_get_data_filtDict(self):
        tmp = [x for x in
               self.img.iter_vector(vector=self.vec, filter={'id':2})]
        self.assertEqual(len(tmp),2)
        self.assertIsNone(tmp[0],None)
        self.assertIsInstance(tmp[1],np.ndarray)

    def test_get_data_filtDictBadKey(self):
        with warnings.catch_warnings(record=True) as w:
            for x in self.img.iter_vector(vector=self.vec, filter={'idbad':2}):
                self.assertIsNone(x)
                self.assertTrue( len(w) == 1 )

    def test_get_data_filtDictBadValue(self):
        for x in self.img.iter_vector(vector=self.vec, filter={'id': 2222}):
            self.assertIsNone(x)

    def test_get_data_filtList(self):
        tmp = [x for x in self.img.iter_vector(vector=self.vec,
                                            filter=[{'id': 2},{'id':1}])]
        self.assertEqual(len(tmp), 2)
        self.assertIsInstance(tmp[0], np.ndarray)
        self.assertIsInstance(tmp[1], np.ndarray)

    def test_get_data_filtListBadValue(self):
        tmp = [x for x in self.img.iter_vector(vector=self.vec,
                                        filter=[{'id': 2222}, {'id': 1}])]
        self.assertEqual(len(tmp), 2)
        self.assertIsInstance(tmp[0], np.ndarray)
        self.assertIsNone(tmp[1])

    def test_get_data_propStr_filtList(self):
        tmp = [x for x in self.img.iter_vector(vector=self.vec,
                        properties='teststr',
                        filter=[{'id': 2},{'id': 1}])]
        self.assertTrue(len(tmp) == 2)
        self.assertIsInstance(tmp[0][0], np.ndarray)
        self.assertTrue(len(tmp[0][1]) == 1)
        self.assertIsInstance(tmp[1][0], np.ndarray)
        self.assertTrue(len(tmp[1][1]) == 1)

    def test_get_data_propList_filtList(self):
        tmp = [x for x in self.img.iter_vector(vector=self.vec,
                          properties=['teststr','testfloat'],
                          filter=[{'id': 2},{'id': 1}])]
        self.assertEqual(len(tmp), 2)
        self.assertIsInstance(tmp[0][0], np.ndarray)
        self.assertTrue(len(tmp[0][1]) == 2)
        self.assertIsInstance(tmp[1][0], np.ndarray)
        self.assertTrue(len(tmp[1][1]) == 2)

class TestDGImage(unittest.TestCase):
    """Testing for :
    DGImage class
    """
    # Run defined method setUp to setup the test environment
    def setUp(self):
        self.test_img = dgsamples.wv2_longmont_1k.ms
        #self.test_img = "../data/imagefiles/053792616010_01/053792616010_01_P001_MUL/14JUN20181517-M2AS-053792616010_01_P001.TIL"
        #self.test_geoimg = "../data/imagefiles/smalldgdelivery/053792616010_01/053792616010_01_P001_MUL/14JUN20181517-M2AS-053792616010_01_P001.TIF"
        #self.test_img_at_sensor_rad = "../data/imagefiles/smalldgdelivery_ms_scaled_to_rad_uW_per_cm2_nm_sr.tif"
        #self.test_img_toa_ref = "../data/imagefiles/smalldgdelivery_ms_scaled_to_toa_ref.tif"
        self.img = geoio.dg.DGImage(self.test_img)
        self.geoimg = geoio.GeoImage(self.test_img)

    def tearDown(self):
        self.img = None
        self.geoimg = None
        # # Remove VRT files from test data dir
        # test_dir = os.path.dirname(os.path.realpath(__file__))
        # search_dir = test_dir+os.sep+os.pardir+os.sep+'data'
        # vrt_files=tt.files.search(search_dir,'*.VRT',
        #                           case_sensitive=False,depth=10)

    def test_DGImage_dg_meta_exists(self):
        self.assertIsInstance(self.img.meta_dg,tt.bunch.OrderedBunch)

    def test_DGImage_dg_meta_abscal(self):
        self.assertTrue(len(self.img.meta_dg_quick.effbandwidth)==8)
        self.assertTrue(len(self.img.meta_dg_quick.abscalfactor)==8)
        self.assertIsInstance(self.img.meta_dg_quick.effbandwidth[0],float)

    def test_DGImage_get_data_as_at_sensor_rad_dtype(self):
        """Output for at sensor radiance should be a numpy array of float32."""
        a = self.img.get_data_as_at_sensor_rad()
        self.assertEqual(a.dtype,np.dtype('float32'))

    def test_DGImage_get_data_as_toa_ref_dtype(self):
        """Output for toa reflectance should be a numpy array of int16 data."""
        a = self.img.get_data_as_toa_ref()
        self.assertEqual(a.dtype,np.dtype('int16'))

    # def test_DGImage_get_data_as_at_sensor_rad_old(self):
    #     """The image loaded into tmp was produced by ENVI."""
    #     a = self.geoimg.get_data()
    #     a = self.img.get_data_as_at_sensor_rad()
    #     tmp = geoio.GeoImage(self.test_img_at_sensor_rad)
    #     b = tmp.get_data()*10  # Convert to W/(m^2*sr*nm)
    #     self.assertIsNone(np.testing.assert_array_almost_equal(a,b,decimal=4))

    def test_DGImage_get_data_as_at_sensor_rad(self):
        """Spot check image against single pixel calculations using the
        following equation:
        L = DN * (ACF/EBW) * (2 - Gain) - Offset
        ACF abscal factor from meta data
        EBW effectiveBandwidth from meta data
        Gain provided by abscal from const
        Offset provided by abscal from const
        """
        a = self.img.get_data_as_at_sensor_rad()
        b = self.img.get_data()
        #b = self.geoimg.get_data()

        # Hard coded from IMD of the test image
        abscalfactor=np.array([0.009295654, 0.01783568,
                               0.01364197, 0.006810718,
                               0.01851735, 0.006063145,
                               0.02050828, 0.009042234])
        effbandwidth=np.array([0.0473, 0.0543, 0.063, 0.0374,
                               0.0574, 0.0393, 0.0989, 0.0996])

        # Dynamic pull from const file that can be updated
        gain=np.array(const.DG_ABSCAL_GAIN['WV02_MULTI'])
        offset=np.array(const.DG_ABSCAL_OFFSET['WV02_MULTI'])

        # Run check for sets of pixel using the simplest per-pixel calc
        for k in range(0,7):
            for i in [100,200,300,400]:
                for j in [100,200,300,400]:
                    pcheck=self.rad_pixel_calc(b,i,j,k,abscalfactor,
                                              effbandwidth,gain,offset)
                    self.assertAlmostEqual(pcheck,a[k,i,j],places=4)

    def rad_pixel_calc(self,data,i,j,k,abscal,effbw,gain,offset):
        # Changed equation Feb. 2016 to match new cal factors from Jan. 2016
        return data[k,i,j]*(abscal[k]/effbw[k])*(gain[k])+offset[k]

    # def test_DGImage_get_data_as_toa_ref_old(self):
    #     """The image loaded into tmp was produced by a running version of my
    #     code with output verified against Python toa_ref code written by
    #     others.
    #     """
    #     a = self.img.get_data_as_toa_ref()
    #     tmp = geoio.GeoImage(self.test_img_toa_ref)
    #     b = tmp.get_data()
    #     self.assertIsNone(np.testing.assert_array_almost_equal(a,b,decimal=6))

    def test_DGImage_get_data_as_toa_ref(self):
        a = self.img.get_data_as_toa_ref()
        b = self.img.get_data()
        #b = self.geoimg.get_data()

        # Hard coded from IMD of the test image
        abscalfactor=np.array([0.009295654, 0.01783568,
                               0.01364197, 0.006810718,
                               0.01851735, 0.006063145,
                               0.02050828, 0.009042234])
        effbandwidth=np.array([0.0473, 0.0543, 0.063, 0.0374,
                               0.0574, 0.0393, 0.0989, 0.0996])
        theta_s=19.20000000000000

        # Dynamic pull from const file that can be updated
        gain=np.array(const.DG_ABSCAL_GAIN['WV02_MULTI'])
        offset=np.array(const.DG_ABSCAL_OFFSET['WV02_MULTI'])

        # Hard coded values for the test image
        e_sun=[1773.81, 2007.27, 1829.62, 1701.85,
               1538.85, 1346.09, 1053.21, 856.599]
        d_es=1.0161415338516235

        # Run check for sets of pixel using the simplest per-pixel calc
        for k in range(0,7):
            for i in [100,200,300,400]:
                for j in [100,200,300,400]:
                    pcheck=self.toa_pixel_calc(b,i,j,k,abscalfactor,
                                              effbandwidth,gain,offset,
                                              e_sun,d_es,theta_s)
                    #print(pcheck)
                    #print(a[k,i,j])
                    self.assertAlmostEqual(pcheck/10000.0,
                                           a[k,i,j]/10000.0,
                                           places=3)

    def toa_pixel_calc(self,data,i,j,k,abscal,effbw,gain,offset,
                       e_sun,d_es,theta_s):
        L=self.rad_pixel_calc(data,i,j,k,abscal,effbw,gain,offset)
        R=(d_es**2*np.pi*L)/(e_sun[k]*np.cos(np.deg2rad(theta_s)))
        return (R*10000).astype('int16')

    def test_DGImage_idl_available(self):
        self.assertTrue(tt.cmd_line.is_tool('idl'))

    def test_DGImage_envi_available(self):
        self.assertTrue(tt.cmd_line.is_tool('envi'))

class TestGeoioFuncs(unittest.TestCase):
    """Testing for :
    Base level functions in Geoio
    """

    def test_request_band_alias_list(self):
        sat_id = 'WV02_MULTI'
        band_alias = ['C','Y']
        out = geoio.dg.get_alias_band_numbers(sat_id, band_alias)
        self.assertEqual(out,[1,4])

    def test_request_band_alias_group(self):
        sat_id = 'WV02_MULTI'
        band_alias = 'VIS'
        out = geoio.dg.get_alias_band_numbers(sat_id, band_alias)
        self.assertEqual(out,[1,2,3,4,5])

    def test_request_band_alias_group2(self):
        sat_id = 'WV02_MULTI'
        band_alias = 'RGBN'
        out = geoio.dg.get_alias_band_numbers(sat_id, band_alias)
        self.assertEqual(out,[5,3,2,7])

    def test_request_band_numbers(self):
        sat_id = 'WV02_MULTI'
        band_alias = [5,4,6]
        out = geoio.dg.get_alias_band_numbers(sat_id, band_alias)
        self.assertEqual(out,[5,4,6])

    def test_request_band_numbers_and_alias(self):
        sat_id = 'WV02_MULTI'
        band_alias = ['C','Y',3]
        out = geoio.dg.get_alias_band_numbers(sat_id, band_alias)
        self.assertEqual(out,[1,4,3])

    def test_request_bad_band_alias_group(self):
        sat_id = 'WV02_MULTI'
        band_alias = 'BADALIAS'
        with self.assertRaises(KeyError):
            geoio.dg.get_alias_band_numbers(sat_id, band_alias)

if __name__ == '__main__':
    unittest.main()