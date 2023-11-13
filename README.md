# imaging_wsclean
This script allows to automatically produce full and low (--dotaper 10'',15'',30'' and --dotaperkpc 25,50,100 kpc) resolutions using WSClean. It also does automatically the source subtraction given the --sourceLLS parameter. The default robust is -0.5 and pixscale is 1.5 [LOFAR setting]
It works for LOFAR, uGMRT and JVLA.
For uGMRT observations, the number of channels is the same as the SP2B.ms; for JVLA observations, the number of channels is the same as the spw (16).
It requires the singularity as it also runs MakeMask.py to improve the cleaning.
To run in the cluster folder

# LOFAR 
python ~/digennaro/scripts/make_cluster_images.py -i CLUSETRNAME --z REDSHIFT --imsize IMSIZE --array LOFAR --dosub --dotaperkpc --dotaper *ms*

# JVLA Lband
python ~/digennaro/scripts/make_cluster_images.py -i CLUSETRNAME --z REDSHIFT --imsize IMSIZE --imsize 4000 --array JVLA --robust 0 --channelsout=16 --pixelscale=1.0 --niter 40000 --dosub --dotaper --dotaperkpc allconf/*ms [allconf folder when we combine multiple array configurations]

# JVLA Sband
python ~/digennaro/scripts/make_cluster_images.py -i CLUSETRNAME --z REDSHIFT --imsize IMSIZE --imsize 4000 --array JVLA --robust 0 --channelsout=16 --pixelscale=0.5 --niter 40000 --dosub --dotaper --dotaperkpc allconf/*ms [allconf folder when we combine multiple array configurations]

# uGMRT band3
python ~/digennaro/scripts/make_cluster_images.py -i CLUSETRNAME --z REDSHIFT --imsize 4000  --array uGMRT --channelsout=N_SP2B --pixelscale=2.0 --dosub --dotaper --dotaperkpc *ms

# uGMRT band4
python ~/digennaro/scripts/make_cluster_images.py -i CLUSETRNAME --z REDSHIFT  --imsize 4000 --array uGMRT --channelsout=N_SP2B --pixelscale=1.0 --dosub --dotaper --dotaperkpc *ms

# uGMRT band5
python ~/digennaro/scripts/make_cluster_images.py -i CLUSETRNAME --z REDSHIFT  --imsize 4000 --array uGMRT --channelsout=N_SP2B --robust=0 --pixelscale=0.5 *.ms
