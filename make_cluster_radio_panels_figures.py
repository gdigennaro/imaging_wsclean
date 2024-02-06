# SCRIPT TO MAKE RADIO IMAGES AT DIFFERENT FREQUENCIES IN SUBPLOTS USING APLPY
# IT CREATES A FIGURE WITH 1xN (WHERE N IS THE NUMBER OF FREQUECIES) OR 2X4 (TOP ROW: FULL VISIBILITIES; BOTTOM ROW: SOURCE SUBTRACTED; COLUMN 1-4: FULL RES, TAPER 25KPC/10arcsec, 50KPC/15arcsec, 100KPC/30arcsec) PANELS
# To run in the folder where the telescope/cluster is
#
# python /export/data/group-brueggen/digennaro/scripts/make_cluster_radio_panels_figures.py -i --DATADIR --z [--RA] [--DEC] --size [--docircle] [--dodiffuse]
#
# G. Di Gennaro
# Sept 2023


import matplotlib.pyplot as plt
import numpy as np
import aplpy
import six
import matplotlib as mpl
import sys, os, glob
import natsort
#import montage_wrapper as montage
import argparse
import itertools
from astropy.visualization import AsymmetricPercentileInterval, ManualInterval
from astropy.visualization import LogStretch
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.utils.data import get_pkg_data_filename
from reproject import reproject_interp
from palettable.cubehelix  import red_16 as cmap

if not sys.warnoptions:
  import warnings
  warnings.simplefilter("ignore")

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

mpl.rcParams['xtick.direction']='in'
mpl.rcParams['ytick.direction']='in'

def fix_aplpy_fits(aplpy_obj, dropaxis=2):
  """This removes the degenerated dimensions in APLpy 2.X...
  The input must be the object returned by aplpy.FITSFigure().
  `dropaxis` is the index where to start dropping the axis (by default it assumes the 3rd,4th place).
  """
  temp_wcs = aplpy_obj._wcs.dropaxis(dropaxis)
  temp_wcs = temp_wcs.dropaxis(dropaxis)
  aplpy_obj._wcs = temp_wcs

def drop2axes(filename, outname):
  hdu = fits.open(filename)[0]
  for kw in "CTYPE", "CRVAL", "CRPIX", "CDELT", "CUNIT":
    for n in 3, 4:
      hdu.header.remove(f"{kw}{n}")
  fits.writeto(outname, hdu.data[0,0], hdu.header, overwrite=True)
 
def computerms(ms,masksup=1.e-7):
  m = ms[np.abs(ms)>masksup]
  rmsold = np.std(m)
  diff = 1.e-1
  cut = 3.
  med = np.median(m)
  for i in range(10):
    ind = np.where(np.abs(m-med)<rmsold*cut)[0]
    rms = np.std(m[ind])
    if np.abs((rms-rmsold)/rmsold)<diff: break
    rmsold = rms
  return rms	


def makeradiofigure(fitsname, z, radec, clustername, IMAGEDIR, Mpcwidth=[0.5,0.5], docircle=True):
  oneradinmpc = cosmo.angular_diameter_distance(z)/(360./(2.*np.pi))
  scalebarlengthdeg    = 1.0/oneradinmpc.value
    
  width = [Mpcwidth[0]/oneradinmpc.value, Mpcwidth[1]/oneradinmpc.value]

  print (clustername, z, radec, width)

  imax = int(len(fitsname))
  irow = int(round(imax/4)) # 4 cols

  if imax <= 5:
    irow, icol = 1, imax
  else:
    icol = int(imax/irow)

  figs = plt.figure( figsize=(icol*4, irow*4))
  for i in range(0, imax):
    print (fitsname[i])

    hdulist = fits.open(fitsname[i])
    rms = computerms(np.ndarray.flatten(hdulist[0].data))  

    # define limits for the image colorbar
    vmin, vmax = rms, 500*rms #100*rms
    
    # define subplot limits
    dx, dy = 0.8/icol, 0.7/irow
    if irow == 1:
      xin, yin = 0.1 + ( (dx+0.028) * (i - (icol * int(i/icol))) ), 0.15 + ( int(i/icol)*(dy+0.05) ) 
    else:
      xin, yin = 0.08 + ( (dx+0.028) * (i - (icol * int(i/icol))) ), 0.13 + ( int(i/icol)*(dy+0.05) )
      
    if i == 0:
      rahide = dechide = False
    elif i in np.arange(0, imax, icol):
      rahide, dechide = True, False
    elif i in np.arange(1, icol):
      rahide, dechide = False, True
    else:
      rahide, dechide = True, True

    print (i, xin, yin, dx, dy)
 

    # APLPY
    fig = aplpy.FITSFigure(fitsname[i], slices=[0,0], figure=figs, subplot=[xin, yin, dx, dy])
    fix_aplpy_fits(fig)
    
    fig.show_colorscale(vmin=vmin, vmax=vmax, stretch='log', cmap=cmap.mpl_colormap, smooth=1)
    
    fig.axis_labels.set_xtext('Right Ascension (J2000)')
    fig.axis_labels.set_ytext('Declination (J2000)')
    fig.axis_labels.set_font(size=16)
    if rahide and dechide:
      fig.axis_labels.hide()
    elif rahide:
      fig.axis_labels.hide_x()
    elif dechide:
      fig.axis_labels.hide_y()

    fig.ticks.set_color('white')
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm')
    fig.tick_labels.set_font(size=13)
    if rahide and dechide:
      fig.tick_labels.hide()
    elif rahide:
      fig.tick_labels.hide_x()
    elif dechide:
      fig.tick_labels.hide_y()

    if False:
      fig.add_grid()
      fig.grid.set_linestyle(':')
      fig.grid.set_color('white')
      fig.grid.set_alpha(0.6)
    
    if radec == [None, None]:
      radec = [float(hdulist[0].header['CRVAL1']), float(hdulist[0].header['CRVAL2'])]
    fig.recenter(radec[0], radec[1], width=width[0], height=width[1])

    fig.set_nan_color('white')

    if docircle:
      fig.show_circles(radec[0],radec[1],0.5*scalebarlengthdeg,edgecolor="white",linestyle=(0,(10,5)),linewidth=1.)
      fig.show_markers(radec[0],radec[1],marker='x',s=30,facecolor='white',linewidth=0.75)   
    
    if i == 2:
      try:
        #name = hdulist[0].header['OBJECT'].strip()
        fig.show_regions(clustername+"_radio_features_im.reg")
      except:
        print ("no region file")
      #  continue

    scale = Mpcwidth[0]/4.
    fig.add_scalebar(scale*scalebarlengthdeg, str(int(scale*1e3))+" kpc",color="white",corner="bottom right",linewidth=1.5,fontsize=10) 

    fig.add_beam(facecolor='crimson',edgecolor='crimson',corner='bottom left',frame=True,alpha=0.7)
    bmaj = round(hdulist[0].header['BMAJ'],4)*3600
    bmin = round(hdulist[0].header['BMIN'],4)*3600
    bpa  = round(hdulist[0].header['BPA'],4)

    telescope = hdulist[0].header['TELESCOP'].strip()
    if telescope == "GMRT": telescope = "uGMRT"
    freq = str(round(hdulist[0].header['CRVAL3']/1.e6))
    
    print ("RESOLUTION [arcsecXarcxsec, deg]:", round(bmaj,1), "x", round(bmin,1), round(bpa))
    print ("NOISE CONTOURS", rms*1e6, "microJy/beam")
    
    # SUBPLOT TITLES
    try:
      fitsname[i].split("sub")[1].split("-MFS-image.fits")[0]
      try:
        taper = fitsname[i].split("TAPER")[1].split("-MFS-image.fits")[0]
        title = freq+" MHz TAPER="+taper+"\nsource subtracted"
      except:
        title = freq+" MHz"+"\nsource subtracted"

    except:
      try:
        taper = fitsname[i].split("TAPER")[1].split("-MFS-image.fits")[0]
        try: 
          fitsname[i].split("kpc")[1]
        except:
          taper = taper+"$''$"
        title = freq+" MHz TAPER="+taper
      except:
        title = freq+" MHz"
    fig.set_title(telescope+" "+title, fontsize=14)

    if True:
      lev_factor = 3.*np.array([-1,1.,2.,4.,8.,16.,32])
      levelsr = np.ndarray.tolist(lev_factor*rms)    
      try:
        fig.show_contour(fitsname[i]+'.contours', slices=[0,0], dimensions=[0,1], levels=levelsr, colors='lightgray', smooth=3, overlap=True,linewidths=1.,alpha=0.5)
      except:  
        drop2axes(fitsname[i], fitsname[i]+'.contours')
        fig.show_contour(fitsname[i]+'.contours', slices=[0,0], dimensions=[0,1], levels=levelsr, colors='lightgray', smooth=3, overlap=True,linewidths=1.,alpha=0.5)

    if docolorbar:
      fig.add_colorbar()
      fig.colorbar.set_location('right')
      fig.colorbar.set_font(size=8) # tick font
      cbar_factor = np.array([1.,2.,4.,8.,16., 32.,64.]) #check your own best ticks
      levelcbar = np.ndarray.tolist(cbar_factor*rms)
      fig.colorbar.set_ticks(levelcbar)
      cbarlabel = [r"$1\rm\sigma_{rms}$", r"$2\rm\sigma_{rms}$", r"$4\rm\sigma_{rms}$", r"$8\rm\sigma_{rms}$", r"$16\rm\sigma_{rms}$", r"$32\rm\sigma_{rms}$", r"$64\rm\sigma_{rms}$"]
      fig.colorbar._colorbar.ax.set_yticklabels(cbarlabel)

  
  #name = hdulist[0].header['OBJECT'].strip()
  figs.suptitle(clustername, fontsize=20)
  figs.tight_layout()

  figs.savefig(imagename+".pdf")
  figs.savefig(imagename+".png") 
  plt.show()
  plt.close()
  return

## MAIN SCRIPT
parser = argparse.ArgumentParser(description='Make radio figures of galaxy cluster(s) for different tapers. We need to have the FITS file in a folder named as the clustername.')

parser.add_argument('-i','--clustername', help='Name of cluster, including band default=image', required=False, type=str)
parser.add_argument('--DATADIR', help='directory where to find the cluster folders',required=True, type=str)
parser.add_argument('--z', help='redshift of cluster',required=True, type=float)
parser.add_argument('--RA', help='RA of cluster in deg', required=False, type=float)
parser.add_argument('--DEC', help='DEC of cluster in deg', required=False, type=float)
parser.add_argument('--size', help='Size of the optical image in Mpc', default=1.0, required=True, type=float)
parser.add_argument('--addcircle', help='in addition to normal imaging, also makes tapered images in uv units', action='store_true')
parser.add_argument('--addcolorbar', help='in addition to normal imaging, also makes tapered images in uv units', action='store_true')
parser.add_argument('--doalltelescopes', help='panels with all frequencies', action='store_true')
parser.add_argument('--dotaperonly', help='only low resolution', action='store_true')
parser.add_argument('--dodiffuse', help='full and low resolution', action='store_true')

args = vars(parser.parse_args())
clustername	= args['clustername']
z					  = args['z']
ra 				  = args['RA']
dec					= args['DEC']
width				= args['size']
DATADIR     = args['DATADIR']+'/*/'+clustername+'/'

imagelist = []
if not args['doalltelescopes']: # single frequency w taper
  imagelist = [sorted(glob.glob(DATADIR+"*_maskROBUST*0*uvmin*-MFS-image.fits"))[0]] +\
                natsort.natsorted(glob.glob(DATADIR+"*_maskROBUST*0*uvmin*TAPER*-MFS-image.fits"))

  if args['dodiffuse']: # single frequency - source subtracted image w tapered
    imagelist   = [sorted(glob.glob(DATADIR+"*_submaskROBUST*0*uvmin*-MFS-image.fits"))[0]] +\
                  natsort.natsorted(glob.glob(DATADIR+"*_masksubROBUST*0*uvmin*TAPER*-MFS-image.fits"))
    
    imagename = args['DATADIR']+"/images/%s_allresolution_diffuse"%clustername
  else:
    imagename = args['DATADIR']+"/images/%s_allresolution"%clustername

else: #for N number of frequencies w/o taper
  imagelist = glob.glob(args['DATADIR']+"/uGMRT/"+clustername[0:5]+"_band3/*_maskROBUST-0.5uvminNone-MFS-image.fits") +\
                glob.glob(args['DATADIR']+"/uGMRT/"+clustername+"_band4/*_maskROBUST-0.5uvminNone-MFS-image.fits") +\
                glob.glob(args['DATADIR']+"/MeerKAT/"+clustername+"/*_maskROBUST-0.5uvminNone-MFS-image.fits")

  if args['dotaperonly']:
  imagelist   = natsort.natsorted(glob.glob(args['DATADIR']+"/uGMRT/"+clustername[0:5]+"_band3/*_maskROBUST*0*uvminNoneTAPER??-MFS-image.fits")) +\
                natsort.natsorted(glob.glob(args['DATADIR']+"/uGMRT/"+clustername+"_band4/*_maskROBUST*0*uvminNoneTAPER??-MFS-image.fits")) +\
                natsort.natsorted(glob.glob(args['DATADIR']+"/MeerKAT/"+clustername[0:5]+"N/*_maskROBUST-0.5uvminNoneTAPER??-MFS-image.fits"))
                   
  imagename = args['DATADIR']+"/images/%s_alltelescope_onlytaper"%clustername

  
  elif args['dodiffuse']: #for N number of frequencies w taper
    imagelist   = glob.glob(args['DATADIR']+"/MeerKAT/"+clustername+"/*_maskROBUST-0.5uvminNone-MFS-image.fits") +\
                  natsort.natsorted(glob.glob(args['DATADIR']+"/MeerKAT/"+clustername+"/*_maskROBUST-0.5uvminNoneTAPER??-MFS-image.fits")) +\
                  glob.glob(args['DATADIR']+"/uGMRT/"+clustername+"_band4/*_maskROBUST-0.5uvminNone-MFS-image.fits") +\
                  natsort.natsorted(glob.glob(args['DATADIR']+"/uGMRT/"+clustername+"_band4/*_maskROBUST-0.5uvminNoneTAPER??-MFS-image.fits")) +\
                  glob.glob(args['DATADIR']+"/uGMRT/"+clustername[0:5]+"_band3/*_maskROBUST-0.5uvminNone-MFS-image.fits") +\
                  natsort.natsorted(glob.glob(args['DATADIR']+"/uGMRT/"+clustername[0:5]+"_band3/*_maskROBUST-0.5uvminNoneTAPER??-MFS-image.fits"))

    imagename = args['DATADIR']+"/images/%s_alltelescope_diffuse"%clustername
  else:
    imagename = args['DATADIR']+"/images/%s_alltelescope"%clustername


#imagelist = allsources + diffuse

print (imagelist)

makeradiofigure(imagelist, z, [ra,dec], clustername, imagename, Mpcwidth=[width,width], docircle=args['docircle'])



sys.exit()
