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
    
if len(filelist) > 0:
   print('===> Identified ', len(filelist), ' Gitm Binary Files.\n')
   testfile = []
   for file in filelist:
       print(file)
       testfile.append(file)
else:
    print("!!! Couldn't find any 3DALL*.bin GITM files !!!!")
    raise Exception()
testfile.sort()
gdata = gitm.GitmBin(testfile[len(filelist) - 1])

PI = np.pi

OrignLons = gdata.attrs['nLon']
OrignLats = gdata.attrs['nLat']
OrignAlts = gdata.attrs['nAlt']

# Pull out actual data
nGCs = 3
index1 = 3
alts = gdata['Altitude' ][0     ,0     ,0:OrignAlts]/1000.0    # Alts in km
lats = gdata['Latitude' ][0     ,nGCs:OrignLats-nGCs,0     ]  # Latitude in radians
lons = gdata['Longitude'][nGCs:OrignLons-nGCs,0     ,0     ]  # Longitude in radians

latsdeg = lats*180.0/PI
lonsdeg = lons*180.0/PI


nLons = len(lons)
nLats = len(lats)
nAlts = len(alts)


AltString = list()

pdfile = 'GlobalMean_Temperature.pdf'

sinlats = np.abs(np.sin(lats))
sumsinlats=np.sum(sinlats)

Temperature = gdata['Temperature'   ][nGCs:OrignLons-nGCs,nGCs:OrignLats-nGCs,:nAlts]  
U           = gdata['V!Dn!N (east)' ][nGCs:OrignLons-nGCs,nGCs:OrignLats-nGCs,:nAlts]
V           = gdata['V!Dn!N (north)'][nGCs:OrignLons-nGCs,nGCs:OrignLats-nGCs,:nAlts]
W           =    gdata['V!Dn!N (up)'][nGCs:OrignLons-nGCs,nGCs:OrignLats-nGCs,:nAlts]

SumLonTemp = np.sum(Temperature, axis=0)/nLons   # Our Zonal Mean
SumLonW    = np.sum(W          , axis=0)/nLons   # Our Zonal Mean

array_shape = (len(lats), len(alts))

#LatWeightedTemp = np.zeros(array_shape,order='F')
GMTemp = np.zeros(nAlts)
GMVelUp = np.zeros(nAlts)

for iAlt in range(nAlts):
    LatWeightedVar = sinlats[0:nLats]*SumLonTemp[0:nLats,iAlt]
    GMTemp[iAlt] = np.sum(LatWeightedVar)/sumsinlats
    LatWeightedVar = sinlats[0:nLats]*SumLonW[0:nLats,iAlt]
    GMVelUp[iAlt] = np.sum(LatWeightedVar)/sumsinlats


# Some plot parameters
axis_text_size = 18
title_text_size = 28
tick_font_size = 13
TitleString = 'Global Man Temperature (K) ' 

fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,11),sharex=True)
line1, = ax.plot(GMTemp, alts)
ax.tick_params(axis='both',labelsize=tick_font_size)    
ax.set_title(TitleString,size = title_text_size)   
ax.set_xlabel('Temperature (K)',size = axis_text_size)
ax.set_ylabel('Altitude (km)',size = axis_text_size)

plt.show()
plt.close()


TitleString = 'Global Mean Velocity (m/s) ' 

fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,11),sharex=True)
line2, = ax.plot(GMVelUp, alts)
ax.tick_params(axis='both',labelsize=tick_font_size)    
ax.set_title(TitleString,size = title_text_size)   
ax.set_xlabel('Velocity (m/s)',size = axis_text_size)
ax.set_ylabel('Altitude (km)',size = axis_text_size)

plt.show()
plt.close()




