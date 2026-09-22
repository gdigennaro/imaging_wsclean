# SCRIPT TO PRODUCE AUTHOMATICALLY SETS OF (TAPERED) RADIO IMAGES WITH AND WITHOUT COMPACT SOURCES FROM A LIST OF CLUSTER
# G. Di Gennaro
#
# May 2019

'''
run inside singularity
'''

import os, glob, argparse
import numpy as np
import sys
from astropy.io import fits
from auxcodes import separator
import time

parser = argparse.ArgumentParser(description='Run extraction and selfcalibration of clusters in LoTSS; you can give either a catalog (FITS format) or the cluster name')
parser.add_argument('-i','--clustername', help='cluster name, if you want to extract a single cluster', default='', required=False, type=str)
parser.add_argument('--z', help='cluster redshift', required=False, type=float)
parser.add_argument('-c','--catalog', help='Catalog to use from which extract clusters', required=False, type=str)
parser.add_argument('--sourceLLS', help='largest linear scale in Mpc to use to create the compact-only image', default=0.4, type=str)
parser.add_argument('--dofirst', help='Use if it is the first time to run or you want to overwrite previous images', action='store_true')


args = vars(parser.parse_args())


if args['catalog'] and args['clustername']:
  print ("Error: either give a single target name or a cluster catalog")
  sys.exit()
elif args['catalog'] and not args['clustername']:
  print ("use catalog:", args['catalog'])
  data = fits.open(args['catalog'])[1].data
  clusterlist = np.array(data['Name']) 
  zlisst      = np.array(data['z']) 
  try:
    LLSlist     = np.array(data['sourceLLS']) 
  except:
    LLSlist     = [np.array(data['R500kpc'][i]/4) for i in range(len(clusterlist))]
elif args['clustername'] and not args['catalog']:
  clusterlist = [args['clustername']]
  zlisst      = [args['z']]
  LLSlist     = [args['sourceLLS']]
else:
  print ("Error: give a single target name or a cluster catalog.")
  sys.exit()

for i, cluster in enumerate(clusterlist):
  start_time = time.time()

  #print ("1",os.getcwd())
  name  = cluster.replace(' ','')
  z     = zlisst[i]
  LLS   = LLSlist[i]
  print ('LLS TO SUBTRACT', LLS)
  print ("")
  print ("CLUSTER:", name, z)

  try:
    os.chdir(name+'/LOFAR/')
    print (os.getcwd())
  except:
    os.chdir(name)
    print (os.getcwd())

  # this deals with the antenna compression of the new facetselfcal version (??)
  if not os.path.exists('./MS/'):
    os.mkdir('./MS/')
    MSs = sorted(glob.glob("*ms*"))
    #for ms in MSs:
    #  os.system("DP3 msin="+ms+" msout=./MS/"+ms+" msout.storagemanager=dysco msout.uvwcompression=false msout.antennacompression=False steps=[]")


  if (args['catalog'] and os.path.exists(name+"_maskROBUST-0.5uvmin80-MFS-image.fits")) and not (args['dofirst']):
    print (name, "imaging done")
    os.chdir('../../')      

  else: #if not os.path.exists(name+"_maskROBUST-0.5uvmin80-MFS-image.fits"):   
    cmd  = "python make_cluster_images.py "
    cmd += "-i "+ name +" "
    if os.path.exists(name+".ds9.reg"):
      cmd += "-b "+ name+".ds9.reg "
    elif os.path.exists(name+"_009-MFS-image.fits"):
      hdu = fits.open(name+"_009-MFS-image.fits")        
      imsize = str(hdu[0].header["NAXIS1"])
      cmd += "--imsize "+ imsize +" "
    else:
      os.system('define an image size')
      sys.exit()
    cmd += "--z "+ str(z) +" "
    cmd += "--array LOFAR "
    cmd += "--sourceLLS "+ str(LLS) +" " #largest linear scale in Mpc
    cmd += "--dosub "
    #cmd += "--dotaper "
    cmd += "--dotaperkpc "
    try:
      cmd += "./*ms*"
    except:
      cmd += "./MS/*ms*"
    print ("")
    print (cmd)
    os.system(cmd)
    os.chdir('../../')
  
  #else:
  #  print (name, "imaging done")
  #  os.chdir('../../')

  timerun = round((time.time() - start_time)/3600,1)
  separator("%s hr"%(str(timerun)))
