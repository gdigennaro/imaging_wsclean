"""
Script to make radio images with the following telescopes:
LOFAR, uGMRT, JVLA.
To be run inside the singularity.
For more information, open README
"""

import matplotlib
matplotlib.use('Agg')
import os, sys
import numpy as np
import glob
from astropy.io import fits
import pyrap.tables as pt
import os.path
import bdsf
import pyregion
import argparse
import pickle
#import aplpy
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def cleanup():
  if not os.path.exists:
    os.mkdir("./singlechannels/")
  os.system('mv *-00*.fits ./singlechannels/.')
  if not os.path.exists("./masks/"):
    os.mkdir("./masks/")
  os.system("mv ./*mask.fits ./masks/.")
  os.system('rm -rf *_subROBUST*.fits') # remove non-masked images
  os.system('rm -rf *_ROBUST*fits') # remove non-masked images
  os.system('rm -rf *-dirty.fits') # remove all dirty images
  return

def getimsize(boxfile, cellsize=1.5):
  """
  find imsize need to image a DS9 boxfile region
  """
  r = pyregion.open(boxfile)

  xs = np.ceil((r[0].coord_list[2])*1.6*3600./cellsize)
  ys = np.ceil((r[0].coord_list[3])*1.6*3600./cellsize)

  imsize = np.ceil(xs) # // Round up decimals to an integer
  if(imsize % 2 == 1):
      imsize = imsize + 1
  return np.int(imsize)


def compute_uvmin(redshift,sourceLLS):
  '''
  sourceLLS in units of Mpc
  '''
  oneradinmpc = cosmo.angular_diameter_distance(redshift)/(360./(2.*np.pi))
  scalebarlengthdeg    = sourceLLS/oneradinmpc.value
    
  return 1./(scalebarlengthdeg*np.pi/180.)
 
def compute_taper(redshift,taperscale):
    '''
    taperscale in units of kpc
    '''
    oneradinmpc = cosmo.angular_diameter_distance(redshift)/(360./(2.*np.pi))
    taper    = 1e-3*taperscale/(oneradinmpc.value)
    
    return taper*3600

def adjustniter_for_taper(taper, niter):
  if taper < 5:
    return np.int(niter)
  if taper >= 5 and taper < 15:
    return np.int(niter/2)
  if taper >= 15:
    return np.int(niter/4)


 
def makeimage(mslist, imageout, pixsize, imsize, channelsout=6, niter=15000, robust=-0.5, minuv=120, uvtaper=None, multiscale=False, predict=True,fitsmask=None, deepmultiscale=False, cluster_redshift=None, column=None):

  # some setup
  username = os.getlogin()
  if username == 'rvweeren':
    wsclean = '/net/lofar1/data1/sweijen/software/LOFAR/latest/wsclean/bin/wsclean'
  else:
    wsclean = 'wsclean'


  if uvtaper != None:
    if uvtaper < 0:
      print ('Not supported uvtaper', uvtaper)
    else:
      imsizein  =  np.int(imsize*(pixsize/(uvtaper/5.)))
      pixsizein = np.int(uvtaper/5.)
      if float(pixsizein) < pixsize: # to deal with rounding issues which cause a 0arcsec pixelsize
        pixsizein = pixsize
        imsizein  = imsize

  else:
    imsizein  = imsize 
    pixsizein = pixsize


  if imsizein < 511: # otherwise images get too small for multiscales
    imsizein = 512


  baselineav = 2.5e3*60000.*2.*np.pi *float(pixsizein)/(24.*60.*60*float(imsizein)) 

  # limit baseline averaging to 10, fixes prob
  if baselineav > 10.0:
    baselineav = 10.

  baselineav = str (baselineav)


  # few simple checks to make sure we have useful data
  msliststring = ' '.join(map(str, mslist))
  os.system('rm -f ' + imageout + '-*.fits')
  imcol = 'CORRECTED_DATA'
  print (mslist[0])
  t = pt.table(mslist[0],readonly=True) # just test for first ms in mslist
  colnames =t.colnames()
  if 'CORRECTED_DATA' not in colnames: # check which column to image
    imcol = 'DATA'
  t.close()
    
  if column != None:
    imcol = column


    
  # build wsclean command      
  cmd = wsclean + ' '
  if minuv:
    cmd += '-minuv-l ' + str(minuv) + ' '
  
  cmd += '-size ' + str(imsizein) + ' ' + str(imsizein) + ' -reorder '
  cmd += '-weight briggs ' + str(robust) + ' '    
  cmd += '-weighting-rank-filter 3 -clean-border 1 '
  cmd += '-mgain 0.8 -fit-beam -data-column ' + imcol + ' '
  cmd += '-join-channels -channels-out ' +str(channelsout) + ' -padding 1.4 '
  #cmd += '-parallel-deconvolution ' + str(np.int(imsizein)/2) + ' ' 
  if multiscale:
    if predict:
      cmd += '-multiscale '+' -multiscale-scales 0,4,8,16 '
    else:
      cmd += '-multiscale '+' -multiscale-scales 0,4,8,16,32,64 '
  
  if fitsmask != None:
    if os.path.isfile(fitsmask): 
      cmd += '-fits-mask '+ fitsmask + ' '
    else:
      print ('fitsmask: ', fitsmask, 'does not exist')
      sys.exit()
  else:
    #cmd += '-auto-mask 1.0 -auto-threshold 0.5 '
    #cmd += '-auto-mask 2.0 -auto-threshold 1.0 '
    cmd += '-auto-mask 3.0 -auto-threshold 1.5 '

  if uvtaper != None:
      cmd += '-taper-gaussian ' +  str(uvtaper) + 'arcsec '
  
  if telescope == "LOFAR":
    cmd += '-fit-spectral-pol 3 '# -beam-shape 6arcsec 6arcsec 0deg '   
    cmd += '-pol I '
    cmd += '-baseline-averaging ' + baselineav + ' '
    cmd += '-no-update-model-required ' 
  elif telescope == "uGMRT":
    cmd += '-pol RR '
    cmd += '-gridder wstacking '
    cmd += '-apply-primary-beam '
  elif telescope == "JVLA":
    cmd += '-apply-primary-beam '
    cmd += '-pol I '
    cmd += '-gridder wstacking '
  elif telescope == 'MeerKAT':
    cmd += '-pol I '
    cmd += '-gridder wstacking '
    
  cmd += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec '

  print ('WSCLEAN: ', cmd + '-niter ' + str(niter) + ' ' + msliststring)
  #logging.info(cmd + '-niter ' + str(niter) + ' ' + msliststring)
  os.system(cmd + '-niter ' + str(niter) + ' ' + msliststring)

  if deepmultiscale:
    # predict first to fill MODEL_DATA so we can continue with clean
    if telescope == "JVLA":
      # check model data are not nan
      modellist = '%s*-00??-model.fits'%(imageout)
      for model in sorted(glob.glob(modellist)): 
        modeldata, hdr = fits.getdata(model, header=True)
        idx = np.where(np.isnan(modeldata))
        modeldata[idx] = 0.
        fits.writeto(model, modeldata, hdr, overwrite=True)

    cmdp = wsclean + ' -size '
    cmdp += str(imsizein) + ' ' + str(imsizein) + ' -channels-out ' + str(channelsout) + ' -padding 1.4 -predict '

    cmdp += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' + msliststring
    print ('PREDICT STEP for continue: ', cmdp)
    os.system(cmdp)
      
    # NOW continue cleaning  
    cmd += '-niter ' + str(niter/15) + ' -multiscale -continue ' + msliststring
    print ('WSCLEAN continue: ', cmd)
    os.system(cmd)

  # REMOVE nagetive model components, these are artifacts (only for Stokes I)
  #if idg:
  #  removenegativefrommodel(sorted(glob.glob(imageout +'-????-I-model*.fits')))  # only Stokes I
  #else:    
  #  removenegativefrommodel(sorted(glob.glob(imageout + '-????-model.fits')))

  if predict:
    if telescope == "JVLA":
      # check model data are not nan
      modellist = '%s*-00??-model.fits'%(imageout)
      for model in sorted(glob.glob(modellist)): 
        modeldata, hdr = fits.getdata(model, header=True)
        idx = np.where(np.isnan(modeldata))
        modeldata[idx] = 0.
        fits.writeto(model, modeldata, hdr, overwrite=True)


    cmd = wsclean + ' -size '
    cmd += str(imsizein) + ' ' + str(imsizein) + ' -channels-out ' + str(channelsout) + ' -padding 1.4 -predict '
    cmd += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' + msliststring
    print ('PREDICT STEP: ', cmd)
    os.system(cmd)
 
 
def subtractcompact(mslist, imageout, pixsize, imsize, minuv, channelsout=6, niter=25000, robust=-0.5, outcolumn='DIFFUSE_SUB'):
  # some setup
  username = os.getlogin()
  if username == 'digennaro':
    makemask = '/net/para10/data1/shimwell/software/killmsddf/new-install/DDFacet/SkyModel/MakeMask.py'
  else:
    makemask = 'MakeMask.py'
  
  makeimage(mslist, imageout +'_compact', pixsize, imsize, channelsout=channelsout, niter=niter, robust=robust, minuv=minuv, multiscale=True, predict=False)

  # make a mask
  imagename  = imageout +'_compact' + '-MFS-image.fits'
  cmdm  = makemask + ' --Th=1.0 --RestoredIm=' + imagename
  print (cmdm)
  os.system(cmdm)
  fitsmask = imagename + '.mask.fits'

  # re-image with mask
  makeimage(mslist, imageout +'_compactmask', pixsize, imsize, channelsout=channelsout, niter=niter, robust=-0.5, minuv=minuv, multiscale=True, predict=True, fitsmask=fitsmask, deepmultiscale=False)

  # now subtract the columns 
  if outcolumn == 'DIFFUSE_SUB':
    for ms in mslist:
      ts  = pt.table(ms, readonly=False)
      colnames = ts.colnames()
      if outcolumn not in colnames:
        desc = ts.getcoldesc('DATA')
        desc['name']=outcolumn
        ts.addcols(desc)
        ts.close() # to write results

      else:
        print (outcolumn, ' already exists')
        ts.close()

    for ms in mslist:
      ts  = pt.table(ms, readonly=False)
      data = ts.getcol('DATA') 
      model = ts.getcol('MODEL_DATA') 
      ts.putcol(outcolumn,data-model)
      ts.close()
  #  
  #THE JVLA DOESN'T ALLOW THE CREATION OF A NEW DATA COLUMN. WE THEREFORE NEED TO COPY THE DATASET IN A NEW ONE, AND OVERWRITE THE CORRECTED_DATA
  elif outcolumn == 'CORRECTED_DATA':
    for ms in mslist:
      #if 'CORRECTED_DATA' in t.colnames():
      print(ms+': Using CORRECTED_DATA-MODEL_DATA')
      os.system("taql 'update " + ms + " set CORRECTED_DATA=CORRECTED_DATA-MODEL_DATA'")
      
      #else:
      #  print('Using DATA-MODEL_DATA') 
      #  os.system("taql 'update " + msnew + " set CORRECTED_DATA=DATA-MODEL_DATA'") 

  return
 
 
 
# some setup
username = os.getlogin()
if username == 'digennaro':
  makemask = '/net/para10/data1/shimwell/software/killmsddf/new-install/DDFacet/SkyModel/MakeMask.py'
else:
  makemask = 'MakeMask.py'
 
 
parser = argparse.ArgumentParser(description='Make images from extraction run. Requires working version of the DR2-pipeline software and WSClean (Oct 2018 or newer)')
parser.add_argument('-b','--boxfile', help='optional boxfile to set imsize automatically', type=str)
#parser.add_argument('--fitsmask', help='fitsmask for deconvolution, if not provided use automasking', type=str)
parser.add_argument('--imsize', help='image size, you can take it from selfcal.log', type=int)
parser.add_argument('-n', '--niter', help='niter, default=35000', default=35000, type=int)
parser.add_argument('--robust', help='Briggs robust paramter, default=-0.5', default=-0.5, type=float)
parser.add_argument('--channelsout', help='channelsout, default=6', default=6, type=int)
parser.add_argument('--minuv', help='inner uv-cut for image in lambda', type=int)
parser.add_argument('--pixelscale', help='pixels size in arcsec, deafult=1.5', default=1.5, type=float)
parser.add_argument('--sourceLLS', help='size in Mpc of diffuse emission for uvcut, deafult=0.4', default=0.4, type=float)
parser.add_argument('--z', help='redshift of cluster, not required if --nodosub is used', default=-1.0, type=float)
parser.add_argument('-i','--imagename', help='imagename, default=image', required=True, type=str)
parser.add_argument('--maskthreshold', help='threshold for MakeMask.py, default=2.0', default=2.0, type=int)
parser.add_argument('--taper', help='Radio image tapers, default=[10,15,30]', default=[10,15,30], type=list)
parser.add_argument('--taperkpc', help='Radio image tapers kpc, default=[25,50,100]', default=[25,50,100], type=list)
#parser.add_argument('--taperkpc', help='Radio image tapers kpc, default=[100]', default=[25,50,100], type=list)
parser.add_argument('--array', help='Interferometer used for the observation (LOFAR/uGMRT/JVLA)', required=True, type=str)
parser.add_argument('--dotaper', help='in addition to normal imaging, also makes tapered images in uv units', action='store_true')
parser.add_argument('--dotaperkpc', help='in addition to normal imaging, also makes tapered images in kpc units', action='store_true')
parser.add_argument('--dosub', help='in addition to normal imaging, also subtract compact sources', action='store_true')
parser.add_argument('ms', nargs='*', help='msfile(s)')

args = vars(parser.parse_args())
 
mslist    = args['ms']
minuv     = args['minuv']  
pixsize   = args['pixelscale']  
niter     = args['niter']  
imageout  = args['imagename']
robust    = args['robust']
tapers    = args['taper']
taperskpc = args['taperkpc']
imsize    = args['imsize']
telescope = args['array']


#if args['boxfile'] == None and args['imsize'] == None:
  #print ('Incomplete input detected, either boxfile or imsize is required')
  #sys.exit()
#if args['boxfile'] != None and args['imsize'] != None:
  #print ('Wrong input detected, both boxfile and imsize are set')
  #sys.exit()

#if args['boxfile'] != None:
  #imsize   = int(getimsize(args['boxfile'], args['pixelscale']))
#if args['imsize'] != None:
  #imsize = int(args['imsize'])


if args['z'] < 0: # if no redshift provided try to find it automatically
  customSimbad = Simbad()
  customSimbad.add_votable_fields('z_value')
  #obj=customSimbad.query_object( args['imagename'] )
  #print obj
  #sys.exit()
  try:
    obj=customSimbad.query_object( args['imagename'] )
  #except: # in this case the cluster name is recognized by Simbad as no error was raised above

    if obj['Z_VALUE'][0] > 0.0:
      print ('Found redshift',  obj['Z_VALUE'][0])
      args['z'] = obj['Z_VALUE'][0]
      os.system('echo ' + str(args['z']) + ' > redshift-used.log')
    else:
      print ('Warning: Cluster known but redshift not found, assuming z=0.2')
      args['z'] = 0.2
      os.system('echo ' + str(args['z']) + ' > redshift-used.log')
  except:
    print ('Cluster name is not known by Simbad')



if args['dosub']:
  #if  telescope == "LOFAR" or telescope == "uGMRT":
  #  columnsub = "DIFFUSE_SUB"

  if telescope == "JVLA":
    # TO BE ADDED: add check that we are using the *ms.sub measurementsets
    columnsub = "CORRECTED_DATA"
  else:
    columnsub = "DIFFUSE_SUB"
  
  if args['z'] < 0:
    print ('You need provide a redshift, none was given')
    sys.exit()

  
  minuv_forsub = compute_uvmin(args['z'], args['sourceLLS'])

  subtractcompact(mslist, imageout, pixsize, imsize, minuv_forsub, channelsout=args['channelsout'], niter=np.int(niter/1.25), robust=robust, outcolumn=columnsub)
  
  #subtractcompact(mslist, imageout, pixsize, imsize, minuv_forsub, channelsout=args['channelsout'], niter=np.int(niter), robust=robust, outcolumn='DIFFUSE_SUB')

if not args['minuv']:
  if telescope == "LOFAR":
    minuv = 80
  else:
    minuv = None

if args['dosub']:
  #  -----------------------------------------------------------------------------
  #  --- make the standard image robust -0.5 image, compact source subtracted ----
  #  -----------------------------------------------------------------------------
  makeimage(mslist, imageout +'_subROBUST'+str(robust)+'uvmin'+str(minuv), pixsize, imsize, channelsout=args['channelsout'], niter=np.int(niter/1.5), robust=robust, minuv=minuv, column=columnsub, predict=False)

  # make a mask
  imagename  = imageout +'_subROBUST'+str(robust)+'uvmin'+str(minuv) + '-MFS-image.fits'
  cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
  print (cmdm)
  os.system(cmdm)
  fitsmask = imagename + '.mask.fits'

  # re-image with mask
  makeimage(mslist, imageout +'_submaskROBUST'+str(robust)+'uvmin'+str(minuv), pixsize, imsize, channelsout=args['channelsout'], niter=np.int(niter/1.5), robust=robust, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, column=columnsub, deepmultiscale=False)    

  if args['dotaper']:
    for t in np.array(tapers):
      #  --------------------------------------------------------
      #  --- make the taper image, compact source subtracted ----
      #  --------------------------------------------------------

      newniter = adjustniter_for_taper(t, niter)
      
      makeimage(mslist, imageout +'_subROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t), pixsize, imsize, channelsout=args['channelsout'], niter=newniter, robust=robust, minuv=minuv, predict=False, column=columnsub, uvtaper=float(t))

      # make a mask
      imagename  = imageout +'_subROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t) + '-MFS-image.fits'
      cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
      print (cmdm)
      os.system(cmdm)
      fitsmask = imagename + '.mask.fits'

      # re-image with mask
      makeimage(mslist, imageout +'_masksubROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t), pixsize, imsize, channelsout=args['channelsout'], niter=newniter, robust=robust, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False, column=columnsub, uvtaper=float(t))

  
  if args['dotaperkpc']:
    if args['z'] >0. :
      for t in np.array(taperskpc):
        #  ---------------------------------
        #  --- make the taper kpc image ----
        #  ---------------------------------
        newniter   = adjustniter_for_taper(compute_taper(args['z'], float(t)), niter)
        uvtaperkpc = compute_taper(args['z'], float(t))

        makeimage(mslist, imageout +'_subROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t)+'kpc', pixsize, imsize, channelsout=args['channelsout'], niter=newniter, robust=robust, minuv=minuv, predict=False, column=columnsub, uvtaper=uvtaperkpc)

        # make a mask
        imagename  = imageout +'_subROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t)+'kpc' + '-MFS-image.fits'
        cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
        print (cmdm)
        os.system(cmdm)
        fitsmask = imagename + '.mask.fits'

        # re-image with mask
        makeimage(mslist, imageout +'_masksubROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t)+'kpc', pixsize, imsize, channelsout=args['channelsout'], niter=newniter, robust=robust, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False, column=columnsub, uvtaper=uvtaperkpc)
    else:
      print ("no redshift for the conversion in kpc")
      #sys.exit()        

if True:
  ##  --------------------------------
  ##  --- make the standard image ----
  ##  --------------------------------  
  makeimage(mslist, imageout +'_ROBUST'+str(robust)+'uvmin'+str(minuv), pixsize, imsize, channelsout=args['channelsout'], niter=niter, robust=robust, minuv=minuv, predict=False)

  # make a mask
  imagename  = imageout +'_ROBUST'+str(robust)+'uvmin'+str(minuv) + '-MFS-image.fits'
  cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
  print (cmdm)
  os.system(cmdm)

  fitsmask = imagename + '.mask.fits'

  # re-image with mask
  makeimage(mslist, imageout +'_maskROBUST'+str(robust)+'uvmin'+str(minuv), pixsize, imsize, channelsout=args['channelsout'], niter=niter, robust=robust, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False)


if False: 
  #  --------------------------------------------------
  #  --- make the high-res image robust -2.0 image ---
  #  --------------------------------------------------
  makeimage(mslist, imageout +'_maskROBUST-2.0'+'uvmin'+str(minuv), pixsize, imsize, channelsout=args['channelsout'], niter=niter, robust=-2.0, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False)

if True: 
  #  --------------------------------------------------
  #  --- make the high-res image robust -1.25 image ---
  #  --------------------------------------------------
  makeimage(mslist, imageout +'_maskROBUST-1.25'+'uvmin'+str(minuv), pixsize, imsize, channelsout=args['channelsout'], niter=niter, robust=-1.25, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False)

 
if args['dotaper']:
  for t in np.array(tapers):
    #  --------------------------------------
    #  --- make the standard taper image ----
    #  --------------------------------------
    newniter = adjustniter_for_taper(t, niter)

    makeimage(mslist, imageout +'_ROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t), pixsize, imsize, channelsout=args['channelsout'], niter=newniter, robust=robust, minuv=minuv, predict=False, uvtaper=float(t))

    # make a mask
    imagename  = imageout +'_ROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t) + '-MFS-image.fits'
    cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
    print (cmdm)
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

    # re-image with mask
    makeimage(mslist, imageout +'_maskROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t), pixsize, imsize, channelsout=args['channelsout'], niter=newniter, robust=robust, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False, uvtaper=float(t))



if args['dotaperkpc']:
  if args['z'] >0. :
    for t in np.array(taperskpc):
      #  ---------------------------------
      #  --- make the taper kpc image ----
      #  ---------------------------------
      newniter = adjustniter_for_taper(compute_taper(args['z'], float(t)), niter)
      uvtaperkpc = compute_taper(args['z'],float(t))

      makeimage(mslist, imageout +'_ROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t)+'kpc', pixsize, imsize, channelsout=args['channelsout'], niter=newniter, robust=robust, minuv=minuv, predict=False, uvtaper=uvtaperkpc)

      # make a mask
      imagename  = imageout +'_ROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t)+'kpc' + '-MFS-image.fits'
      cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
      print (cmdm)
      os.system(cmdm)
      fitsmask = imagename + '.mask.fits'

      # re-image with mask
      makeimage(mslist, imageout +'_maskROBUST'+str(robust)+'uvmin'+str(minuv)+'TAPER'+str(t)+'kpc', pixsize, imsize, channelsout=args['channelsout'], niter=newniter, robust=robust, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False, uvtaper=uvtaperkpc)
  else:
    print ("no redshift for the conversion in kpc")
    #sys.exit()

cleanup()

