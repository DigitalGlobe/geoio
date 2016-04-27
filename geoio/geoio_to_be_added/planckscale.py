#!/usr/bin/env python

from osgeo import gdal
from osgeo import gdalconst
import shutil
import os
from optparse import OptionParser
import numpy as np
import ipdb
from IPython import embed



def planck(wav,T):
  """Return Plank black body curve.  Input wavelength in nanometers, temp in Kelvin."""
  h = 6.626e-34
  c = 3.0e+8
  k = 1.38e-23
  
  wav = np.asarray(wav,dtype=float)*10**-9
  
  return (2*h*c**2)/(wav**5*(np.exp(h*c/k/T/wav)-1))

def find_scale_factors(): 
  wav = np.arange(350,2450+1)
  T = 5778.0  # Sun temp in K 
  
  bbcurve = planck(wav,T)

  # Read wv3 normalized spectral curves
  wv3_curves = np.loadtxt(
            '/nlongbot/rss/users/nwl/swir_vnir_3p7_samples/wv3_norm_spectral_nm.txt',
            skiprows=18).T

  # Multiply with wv3 normalized spectral curves
  bb_mean = []
  bb_max = []
  for i in range(16): bb_mean.append((bbcurve/bbcurve.max()*wv3_curves[i+1,:]).sum())
  for i in range(16): bb_max.append((bbcurve/bbcurve.max()*wv3_curves[i+1,:]).max())

  return bb_max

if __name__=="__main__":

  # parse command line
  parser = OptionParser()
  parser.add_option('--inputfile',help='16 band file to copy and scale',default=None)
  (args,args_extra) = parser.parse_args()
  
  # Figure new file name and copy
  fin = args.inputfile
  (root,ext) = os.path.splitext(fin)
  fname = root+"_scaled"+ext 
  print "Copying file from:"
  print fin
  print "to:"
  print fname
  shutil.copyfile(fin,fname)
  #fname = "/nlongbot/rss/users/nwl/swir_vnir_3p7_samples/goldfield/goldfield_1p6_from3p7_stack_scaled.tif"
  
  # Open the vnir/swir stack
  f = gdal.Open(fname,gdalconst.GA_Update)
  
  # Set (hard code for now) the scaling factors for each band
  bottom = [1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000]
  top = [5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000]
  mean = [2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000]

  scale = find_scale_factors()

  # Scale each band and save it back to disk
  for i in range(f.RasterCount):
    print 'Starting band '+str(i+1)+' update...'
    b = f.GetRasterBand(i+1)
   
    ba = np.ma.masked_less_equal(b.ReadAsArray(),0.0)
    bam = np.ma.asarray(ba,dtype=float)
    
    ttt = scale[i]/max(scale) * 10000
    bbb = ttt-(bam.max()-bam.min()) 
    print ttt
    print bbb

    bam = ((bam-bam.min())*(ttt-bbb))/(bam.max()-bam.min())+bbb
    bam = np.ma.asarray(bam,dtype='UInt16')
  
    b.WriteArray(bam)
    b = None
  
  # Close the image file to force flush to disk
  f = None
  
