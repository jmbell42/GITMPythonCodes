# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import spacepy
import gitm 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import cmocean
import glob, os
from cartopy import config
import cartopy.crs as ccrs
from matplotlib.backends.backend_pdf import PdfPages
# Pull in the rc from matplotlib to make pretty fonts (LaTeX Fonts)
from matplotlib import rc
from matplotlib import ticker
import matplotlib as mpl
import pdb

# First, we find all 3DALL binary files in the directory and then
# we simply access the first one.

# Set some global parameters for Text
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
mpl.rcParams["lines.linewidth"] = 1.5
mpl.rcParams["contour.negative_linestyle"]='dashed'


# Turn off Interactive MODE (Suppresses the Figures when writing files)
plt.ioff()

# ACCESS GITM DATA FILES
# FILE 1
print('===> Searching Working Directory for GITM 3DALL Files. \n')
filelist = []  # initialize a filelist array
for file in glob.glob("3DALL*.bin"):
    filelist.append(file)
    
filelist.sort()
if len(filelist) > 0:
   print('===> Identified ', len(filelist), ' Gitm Binary Files.\n')
   testfile = []
   for file in filelist:
       print(file)
       testfile.append(file)
else:
    print("!!! Couldn't find any 3DALL*.bin GITM files !!!!")
    raise Exception()

gdata = gitm.GitmBin(testfile[len(filelist)-1])

PI = np.pi

OrignLons = gdata.attrs['nLon']
OrignLats = gdata.attrs['nLat']
OrignAlts = gdata.attrs['nAlt']

# Pull out actual data
nGCs = 2
alts = gdata['Altitude' ][0     ,0     ,0:OrignAlts]/1000.0      # Alts in km
lats = gdata['Latitude' ][0     ,0:OrignLats,0     ]*(180.0/PI)  # Lats in deg
lons = gdata['Longitude'][0:OrignLons,0     ,0     ]*(180.0/PI)  # Lons in deg

nLons = len(lons)
nLats = len(lats)
nAlts = len(alts)

for iAlt in range(nAlts-1):
    print('iAlt = ', iAlt, '. Altitude = ', alts[iAlt])
    print('iAlt = ', iAlt, '. dZ     e = ', alts[iAlt+1] - alts[iAlt])
