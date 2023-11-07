# SCRIPT TO PRODUCE AUTHOMATICALLY SETS OF (TAPERED) RADIO IMAGES WITH AND WITHOUT COMPACT SOURCES FROM A LIST OF CLUSTER
# G. Di Gennaro
#
# May 2019

'''
run on singularity:
singularity run --pid --writable-tmpfs --containall --cleanenv -B /export/,/hsopt/ /export/data2/AG_Brueggen/digennaro/software/pill-latest.simg
'''

import os, glob
import numpy as np
import sys
from astropy.io import fits

DATADIR = "/export/data/group-brueggen/digennaro/HighRedshiftClusters/MaDCoWS/LoTSS/"
data = fits.open(DATADIR + "madcows_lofar-dr2.fits")[1].data #np.loadtxt(DATADIR+"madcows_lofar-dr2.txt", dtype="str")

for i, clustername in enumerate(data['clustername']):
  print (clustername)

  z = data['photoz'][i]
  MS = DATADIR + clustername +"/*ms*"

  #if not        
  if not os.path.exists(DATADIR+clustername+"/"+clustername+"_maskROBUST-0.5uvmin80-MFS-image.fits"):   
    if os.path.exists(DATADIR+clustername+"/"+clustername+"_image_9-MFS-image.fits") or os.path.exists(DATADIR+clustername+"/"+clustername+"_image_009-MFS-image.fits"):
      
      fitsimage = glob.glob(DATADIR+clustername+"/"+clustername+"_image_*9-MFS-image.fits")[0]
      hdu = fits.open(fitsimage)        
      imsize = str(hdu[0].header["NAXIS1"])

      cmd  = "python make_cluster_images.py "
      cmd += "-i "+ clustername +" "
      cmd += "--imsize "+ imsize +" "
      cmd += "--z "+ str(z) +" "
      cmd += "--array LOFAR "
      cmd += "--sourceLLS 0.25 "
      cmd += "--dosub "
      #cmd += "--dotaper "
      cmd += "--dotaperkpc "
      cmd += MS

      print ("")
      print (cmd)
      #os.system(cmd)
      #os.system("mv "+clustername+"*fits* "+ clustername +"/.")
      print ("")   
    
    else:
      print ("no imsize available")


  else:
    print (clustername, "done")

