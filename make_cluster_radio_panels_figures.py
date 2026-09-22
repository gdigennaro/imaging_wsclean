# SCRIPT TO PRODUCE AUTHOMATICALLY SETS OF (TAPERED) RADIO IMAGES WITH AND WITHOUT COMPACT SOURCES FROM A LIST OF CLUSTER
# G. Di Gennaro
#
# May 2019
# Modified to support multiple subplots within a single figure (analogous to make_cluster_radio_panels_figures.py)

import os, glob, argparse
import numpy as np
import sys
from astropy.io import fits
import time
import math

import matplotlib.pyplot as plt
import numpy as np
import aplpy
import six
import matplotlib as mpl
import sys, os, glob
import natsort
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
from matplotlib import cm as color_map


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
  fits.writeto(outname, hdu.data[:,:], hdu.header, overwrite=True)

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

def radius(M500, z):
  rhoc = cosmo.critical_density(z)
  M500 *= 1.e14*(2.e33) #M500 in g
  R500 = ( (3.*M500)/(4.*np.pi* (500*rhoc.value)) )**(1./3.) #R500 in cm
  R500Mpc = R500/3.08e24 #R500 in Mpc
  R500deg = round((R500Mpc*1e3 / cosmo.kpc_proper_per_arcmin(z).value / 60.),5) #R500 in deg
  return (R500Mpc, R500deg)


def makeradiofigure(fitsnames, z, radec, M, clustername, scale, f, outname, sizewidth, regions,
                    dosizedeg=False, docircle=True, docolorbar=True, set_panelstitle=False,
                    set_paneltitle=False, nocontours=False, nolabels=False, notick=False, doscale=True, dobeam=True,
                    addextracontours=False, mapcolor='red', dovertical=False):

  # Accept either a single filename (str) or a list
  if isinstance(fitsnames, str):
    fitsnames = [fitsnames]

  oneradinmpc = cosmo.angular_diameter_distance(z) / (360. / (2. * np.pi))
  scalebarlengthdeg = 1.0 / oneradinmpc.value

  if dosizedeg:
    width = [sizewidth[0], sizewidth[1]]
  else:
    width = [sizewidth[0] / oneradinmpc.value, sizewidth[1] / oneradinmpc.value]

  # -------------------------------------------------------------------------
  # Grid layout (mirrors make_cluster_radio_panels_figures.py)
  # -------------------------------------------------------------------------
  imax = len(fitsnames)
  singleplot = (imax == 1)

  if imax <= 5:
    irow, icol = 1, imax
  else:
    irow = int(round(imax / 4))   # 4 columns
    icol = int(imax / irow)

  print(clustername, z, radec, width)

  #figs = plt.figure(figsize=(icol * 6, irow * 5.5))
  if dovertical:
    figs = plt.figure( figsize=(irow*5.5, icol*3.5))
    print("Grid: %d row(s) x %d col(s)" % (icol, irow))  
    dx, dy = 0.8/icol, 0.7/irow
  
  else:
    figs = plt.figure( figsize=(icol*4.4, irow*3.5))
    print("Grid: %d row(s) x %d col(s)" % (irow, icol))
    dx, dy = 0.8/icol, 0.7/irow


  for i, fitsname in enumerate(fitsnames):
    print(fitsname)

    hdulist = fits.open(fitsname)
    rms = computerms(np.ndarray.flatten(hdulist[0].data))

    # Colorbar limits
    if f == [None, None]:
      vmin, vmax = rms, 500 * rms
    else:
      vmin, vmax = f[0] * rms, f[1] * rms

    # ------------------------------------------------------------------
    # Subplot position (same arithmetic as make_cluster_radio_panels_figures.py)
    # ------------------------------------------------------------------

    if irow == 1:
      xin = 0.1  + ((dx + 0.028) * (i - (icol * int(i / icol))))
      yin = 0.15 + (int(i / icol) * (dy + 0.05))
      if icol == 2:
        xin = 0.15 + ((dx + 0.028) * (i - (icol * int(i / icol))))
    else:
      xin = 0.08 + ((dx + 0.028) * (i - (icol * int(i / icol))))
      yin = 0.13 + (int(i / icol) * (dy + 0.05))

    # Axis-label / tick visibility
    col_idx = i % icol
    row_idx = i // icol
    if i == 0:
      rahide = dechide = False
    elif col_idx == 0:          # leftmost column (not first panel)
      rahide, dechide = True, False
    elif row_idx == 0:          # top row (not first panel)
      rahide, dechide = False, True
    else:
      rahide, dechide = True, True

    print(i, xin, yin, dx, dy)

    # ------------------------------------------------------------------
    # APLpy figure
    # ------------------------------------------------------------------
    #try:
    #  drop2axes(fitsname,fitsname)
    #except:
    #  pass
    if singleplot:
      fig = aplpy.FITSFigure(fitsname, slices=[0,0], figure=figs, smooth=1)
    else:
      if dovertical:
        fig = aplpy.FITSFigure(fitsname, slices=[0,0], figure=figs,
                             subplot=[yin, xin, dy, dx])
      else:
        fig = aplpy.FITSFigure(fitsname, slices=[0,0], figure=figs,
                             subplot=[xin, yin, dx, dy])


    try:
      fix_aplpy_fits(fig)
    except:
      pass

    if mapcolor == 'red':
      from palettable.cubehelix  import red_16 as cmap
    elif mapcolor == 'mycubehelix':
      from palettable.cubehelix  import cubehelix2_16 as cmap
    elif mapcolor == 'cubehelix':
      from palettable.cubehelix  import cubehelix_16 as cmap      
    elif mapcolor == 'nicebw':
      color_map_name = 'bone_r'
    elif mapcolor == 'bw':
      from palettable.cmocean.sequential import Gray_8_r as cmap

    if mapcolor != 'nicebw': color_map_name = cmap.mpl_colormap
    contourscolor  = 'lightgray'

    my_map = color_map.get_cmap(color_map_name)
    colors = my_map(np.linspace(0, 1, 3, endpoint=True))

    fig.show_colorscale(vmin=vmin, vmax=vmax, stretch=scale,
                        cmap=color_map_name, smooth=3)

    # Axis labels
    fig.axis_labels.set_xtext('Right Ascension (J2000)')
    fig.axis_labels.set_ytext('Declination (J2000)')
    fig.axis_labels.set_font(size=13)
    if rahide and dechide:
      fig.axis_labels.hide()
    elif rahide:
      fig.axis_labels.hide_x()
    elif dechide:
      fig.axis_labels.hide_y()
    if nolabels:
      fig.axis_labels.hide()

    # Tick labels
    fig.ticks.set_color(contourscolor)
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm')
    fig.tick_labels.set_font(size=12)
    if rahide and dechide:
      fig.tick_labels.hide()
    elif rahide:
      fig.tick_labels.hide_x()
    elif dechide:
      fig.tick_labels.hide_y()
    if notick:
      fig.tick_labels.hide()
      fig.ticks.hide()

    # Recentre
    if radec == [None, None]:
      radec = [float(hdulist[0].header['CRVAL1']),
               float(hdulist[0].header['CRVAL2'])]
      print(360 + radec[0], radec[1])
    fig.recenter(radec[0], radec[1], width=width[0], height=width[1])

    fig.set_nan_color('white')

    # Region files
    if regions:
      for regfile in regions:
        fig.show_regions(regfile)

    # Scale bar
    if doscale:
      scalebar = 0.5  # 500 kpc
      fig.add_scalebar(scalebar * scalebarlengthdeg, "500 kpc",
                       color=contourscolor, corner="bottom",
                       linewidth=1.5, fontsize=10)

    # Beam
    if dobeam:
      fig.add_beam(facecolor=colors[1], edgecolor=colors[1],
                   corner='bottom left', frame=True)
      bmaj = round(hdulist[0].header['BMAJ'] * 3600, 1)
      bmin = round(hdulist[0].header['BMIN'] * 3600, 1)
      bpa  = round(hdulist[0].header['BPA'])
      print("RESOLUTION [arcsecXarcsec, deg]:", bmaj, "x", bmin, bpa)

    if set_paneltitle:
      title = hdulist[0].header['TELESCOP'].strip()
      if title == "GMRT": title = "uGMRT"
   
      try:
        freq = str(round(hdulist[0].header['CRVAL3'] / 1.e6))
      except:
        freq = str(round(hdulist[0].header['FREQ']   / 1.e6))
      title += " " + freq + " MHz"

      if dobeam:
        # Panel title
        if 'sub' in str(fitsname):
          title += r" $\Theta=" + str(round(bmaj,1)) + r"''\times" + str(round(bmin,1)) + "''$\nsource subtracted"
        elif 'compact' in str(fitsname):
          if   'image'    in str(fitsname): suffix = "compact only"
          elif 'model'    in str(fitsname): suffix = "compact only - model"
          elif 'residual' in str(fitsname): suffix = "compact only - residual"
          #else:                             suffix = "compact"
          title += r" $\Theta=" + str(round(bmaj,1)) + r"''\times" + str(round(bmin,1)) + "''$\n" + suffix
        else:
          title += r" $\Theta=" + str(round(bmaj,1)) + r"''\times" + str(round(bmin,1)) + "''$"

      fig.set_title(title, fontsize=11)

    # Contours
    if not nocontours:
      lowrescontours = fitsnames[i] #glob.glob('./%s/LOFAR/*_masksubROBUST-0.5uvmin80TAPER100kpc-MFS-image.fits' % clustername)

      if lowrescontours:
        #lowrescontours = lowrescontours[0]
        hdulistlow = fits.open(lowrescontours)
        rmslow = computerms(np.ndarray.flatten(hdulistlow[0].data))
        print("NOISE CONTOURS", rmslow * 1e6, "microJy/beam")
        lev_factor = 2.5 * np.array([1, 2., 4.])
        levelsr    = np.ndarray.tolist(lev_factor * rmslow)
        try:
          fig.show_contour(lowrescontours, slices=[0,0], dimensions=[0,1],
                           levels=levelsr, colors=contourscolor, smooth=3,
                           overlap=True, linewidths=1, alpha=0.8)
        except:
          try:
            fig.show_contour(lowrescontours + '.contours', slices=[0,0], dimensions=[0,1],
                             levels=levelsr, colors=contourscolor, smooth=3,
                             overlap=True, linewidths=1, alpha=0.8)
          except:
            drop2axes(lowrescontours, lowrescontours + '.contours')
            fig.show_contour(lowrescontours + '.contours', slices=[0,0], dimensions=[0,1],
                             levels=levelsr, colors=contourscolor, smooth=3,
                             overlap=True, linewidths=1, alpha=0.8)
      else:
        print("WARNING: no low-res contour file found; skipping contours for panel %d" % i)

      if addextracontours:
        SZcontours = glob.glob('./%s/SZ/*fits' % clustername)[0]
        print(SZcontours)
        hdulistSZ = fits.open(SZcontours)
        rmsSZ = computerms(np.ndarray.flatten(hdulistSZ[0].data))
        lev_factorSZ = 2. * np.array([1, 2., 4., 6., 8., 10., 16., 32])
        levelsrSZ    = np.ndarray.tolist(lev_factorSZ * rmsSZ)
        fig.show_contour(SZcontours, slices=[0,0], dimensions=[0,1],
                         levels=levelsrSZ, linestyles='--', colors='cyan',
                         smooth=1, overlap=True, linewidths=0.7, alpha=0.65)

    # Colorbar
    if docolorbar:
      fig.add_colorbar()
      fig.colorbar.set_location('right')
      fig.colorbar.set_font(size=8)
      cbar_factor = np.array([1., 2., 4., 8., 16., 32., 64.])
      levelcbar   = np.ndarray.tolist(cbar_factor * rms)
      fig.colorbar.set_ticks(levelcbar)
      cbarlabel = [r"$1\rm\sigma_{rms}$", r"$2\rm\sigma_{rms}$",
                   r"$4\rm\sigma_{rms}$", r"$8\rm\sigma_{rms}$",
                   r"$16\rm\sigma_{rms}$", r"$32\rm\sigma_{rms}$",
                   r"$64\rm\sigma_{rms}$"]
      fig.colorbar._colorbar.ax.set_yticklabels(cbarlabel)

    # R500 circles
    if docircle:
      R500Mpc, R500deg = radius(M, z)
      fig.show_circles(radec[0], radec[1], 0.5 * R500deg,
                       edgecolor="white", linestyle=(0,(8,10,1,10)), linewidth=1)
      fig.add_label(radec[0], radec[1] + (0.5 * R500deg * 1.1),
                    r'0.5$R_{500}$', color='w', size=9, ha='center')
      fig.show_circles(radec[0], radec[1], R500deg,
                       edgecolor="white", linestyle=(0,(10,5)), linewidth=0.75)
      fig.add_label(radec[0], radec[1] + (R500deg * 1.1),
                    r'$R_{500}$', color='w', size=9, ha='center')
      fig.show_markers(radec[0], radec[1], marker='x', s=30,
                       facecolor='white', linewidth=0.75)

  # Overall title
  if set_panelstitle:
    figs.suptitle(clustername +
                  r' [$z=%s$, $M_{500}=%s\times10^{14}~{\rm M_\odot}$]' % (str(z), str(round(M, 1))),
                  fontsize=16)

  figs.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

  if not os.path.exists('./images/'):
    os.mkdir('./images/')
  output = "./images/%s" % outname

  figs.savefig(output + ".pdf", bbox_inches='tight')
  figs.savefig(output + ".png", bbox_inches='tight')
  plt.show()
  plt.close()
  return


# =============================================================================
# CLI
# =============================================================================
parser = argparse.ArgumentParser(
    description='Run extraction and selfcalibration of clusters in LoTSS; '
                'you can give either a catalog (FITS format) or the cluster name')
parser.add_argument('-i', '--clustername', help='cluster name, if you want to extract a single cluster',
                    default='', required=False, type=str)
parser.add_argument('--z',   help='cluster redshift',       required=False, type=float)
parser.add_argument('--RA',  help='cluster RA (in deg)',    required=False, type=float)
parser.add_argument('--DEC', help='cluster DEC (in deg)',   required=False, type=float)
parser.add_argument('--M',   help='cluster mass',           required=False, type=float)
parser.add_argument('-c', '--catalog', help='Catalog to use from which extract clusters',
                    required=False, type=str)
parser.add_argument('--size',  nargs=2, help='Size of the image (Mpc or deg)',
                    required=False, type=float)
parser.add_argument('--dosizedeg', help='if size is in degrees', default=False, action='store_true')
parser.add_argument('--range', nargs=2, help='factor to multiply to rms for vmin/vmax',
                    default=[0.5, 100], required=False, type=float)
parser.add_argument('--mapstretch', help='map color stretch', default='log', type=str)
parser.add_argument('--mapcolor',   help='map color for palettable', required=False, type=str)
parser.add_argument('--dofirst', help='Use if it is the first time to run or to overwrite previous images',
                    action='store_true')
parser.add_argument('--addregion', nargs='+', help='if you want to add a region, use the name', required=False)
parser.add_argument('--nocontours',   help='if you want to suppress contours',  action='store_true')
parser.add_argument('--addcolorbar',  help='add colorbar to each panel',        action='store_true')
parser.add_argument('--addcircle',    help='add R500 circles',                  action='store_true')
parser.add_argument('--addclustertitle', help='add title with cluster name',    action='store_true', default=False)
parser.add_argument('--addpaneltitle', help='add title for each pane',    action='store_true', default=False)
parser.add_argument('--noticks',  help='suppress ticks',  action='store_true')
parser.add_argument('--nolabels', help='suppress axis labels', action='store_true')
parser.add_argument('--doscale',  help='add kpc scale bar', action='store_true')
parser.add_argument('--dobeam',   help='add beam ellipse',  action='store_true')
parser.add_argument('--dovertical',   help='if to make horizontal panels',  action='store_true')
parser.add_argument('--fitsimage', nargs='+', help='FITS image(s)', required=False)
parser.add_argument('-o', '--outname', help='output filename (without extension)', required=False, type=str)

args = vars(parser.parse_args())

if args['catalog'] and args['clustername']:
  print("Error: either give a single target name or a cluster catalog")
  sys.exit()
elif args['catalog'] and not args['clustername']:
  print("use catalog:", args['catalog'])
  data        = fits.open(args['catalog'])[1].data
  clusterlist = np.array(data['Name'])
  zlist       = np.array(data['z'])
  ralist      = np.array(data['RAJ2000'])
  declist     = np.array(data['DEJ2000'])
  Mlist       = np.array(data['M500'])
elif args['clustername'] and not args['catalog']:
  clusterlist = [args['clustername']]
  zlist       = [args['z']]
  ralist      = [args['RA']]
  declist     = [args['DEC']]
  Mlist       = [args['M']]
else:
  print("Error: give a single target name or a cluster catalog.")
  sys.exit()

for i, cluster in enumerate(clusterlist):

  name = cluster.replace(' ', '')
  z    = zlist[i]
  ra   = ralist[i]
  dec  = declist[i]
  if args['M']:
    M = Mlist[i]; docircle = True
  else:
    M = None; docircle = False

  print("")
  print("CLUSTER:", name, z, ra, dec, M)

  if args['size']:
    size = args['size']
  else:
    R500Mpc = radius(M, z)[0]
    size = [2.5 * R500Mpc, 2.5 * R500Mpc]
    print(R500Mpc, size)

  if (args['catalog']) and not (args['dofirst']):
    print("image already exists")
  else:
    # Build the image list: --fitsimage accepts 1..N files
    if args['fitsimage']:
      imagenames = args['fitsimage']   # already a list thanks to nargs='+'
    else:
      imagenames = ["./%s/LOFAR/%s_maskROBUST-0.5uvmin80-MFS-image.fits" % (name, name)]

    outname = args['outname'] if args['outname'] else name + '_panels'

    makeradiofigure(
      imagenames, z, [float(ra), float(dec)], M, name,
      args['mapstretch'], [args['range'][0], args['range'][1]], outname,
      [size[0], size[1]], regions=args['addregion'],
      dosizedeg=args['dosizedeg'],
      docircle=args['addcircle'],
      docolorbar=args['addcolorbar'],
      set_panelstitle=args['addclustertitle'],
      set_paneltitle=args['addpaneltitle'],
      nocontours=args['nocontours'],
      nolabels=args['nolabels'],
      notick=args['noticks'],
      doscale=args['doscale'],
      dobeam=args['dobeam'],
      addextracontours=False,
      mapcolor=args['mapcolor'],
      dovertical=args['dovertical']
    )
