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
import matplotlib as mpl
import pdb 

# First, we find all 3DALL binary files in the directory and then
# we simply access the first one.

#------------------------------------------------------
# Set some global parameters for Text
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
mpl.rcParams["lines.linewidth"] = 1.5
mpl.rcParams["contour.negative_linestyle"]='dashed'
# Turn off Interactive MODE (Suppresses the Figures when writing files)
plt.ioff()


#-----------------------------------------------------
# ACCESS GITM DATA FILES
# FILE 1
print('===> Searching Working Directory for GITM 3DALL Files. \n')
filelist = []  # initialize a filelist array
for file in glob.glob("*.bin"):
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
nGCs = 0
alts = gdata['Altitude' ][0     ,0     ,0:OrignAlts]/1000.0      # Alts in km
lats = gdata['Latitude' ][0     ,3:OrignLats-3,0     ]*(180.0/PI)  # Lats in deg
lons = gdata['Longitude'][3:OrignLons-3,0     ,0     ]*(180.0/PI)  # Lons in deg

nLons = len(lons)
nLats = len(lats)
nAlts = len(alts)

# These are the indices of the Altitudes that we want to pull out.
LonString = list()

#-----------------------------------------------------
# BEGIN OUTPUT OF PDF FILES
#-----------------------------------------------------
axis_text_size = 14
title_text_size = 28
tick_font_size = 13
subpanel_text_size = 14

# Note the .sum(axis=0) will sum over all longitudes(first index)
# we divide by nLons to make it a zonal mean
Temperature = gdata['Temperature'   ][3:OrignLons-3,3:OrignLats-3,:nAlts].sum(axis=0)/nLons
U           = gdata['V!Dn!N (east)' ][3:OrignLons-3,3:OrignLats-3,:nAlts].sum(axis=0)/nLons
V           = gdata['V!Dn!N (north)' ][3:OrignLons-3,3:OrignLats-3,:nAlts].sum(axis=0)/nLons

nN2         = gdata['N2'            ][3:OrignLons-3,3:OrignLats-3,:nAlts].sum(axis=0)/nLons  
nCH4        = gdata['CH4'           ][3:OrignLons-3,3:OrignLats-3,:nAlts].sum(axis=0)/nLons  
nH2         = gdata['N2'            ][3:OrignLons-3,3:OrignLats-3,:nAlts].sum(axis=0)/nLons 
NT = nN2 + nCH4 + nH2
 


#pdb.set_trace()
AltString = list()


pdfile = 'CheckingTemperature_LBC.pdf'
with PdfPages(pdfile) as pdf:
    #for i in AltitudeSlices:
    for i in range(nAlts):
        altpick = alts[i]
        newstring = str(altpick)
        splitstring = newstring.split('.')
        AltString.append(splitstring[0])

    AltSliceIndex = 0
    #for iAlt in range(nAlts):
    for iAlt in range(10):
        print('Accessing iAlt %i',iAlt)
        iAlt = min(iAlt,nAlts-1)
        TitleString = 'Altitude = ' + \
                      AltString[iAlt] + ' km'
        # Figsize = (width, height) in inches
        # subplots (nrows, ncolumns, figure_settings)
        fig, axs = plt.subplots(nrows=2,ncols=2,figsize=(8.5,5.0),sharex=True)

        axs[0,0].annotate(TitleString,xy=(0.125,0.935), \
                    xycoords='figure fraction',\
                    fontsize=title_text_size)
        plotarray = Temperature[:nLats,iAlt]
        plotarray_north = np.ma.masked_where(lats <= 0.0, plotarray)
        lats_north      = np.ma.masked_where(lats <= 0.0, lats)

        plotarray_south = np.ma.masked_where(lats >= 0.0, plotarray)
        lats_south      = np.ma.masked_where(lats >= 0.0, lats)

        #axs[0,0].plot(           lats,plotarray      ,color='black', linestyle = 'solid')
        axs[0,0].plot(    lats_north ,plotarray_north,color='blue' ,linestyle='dashed')
        axs[0,0].plot(abs(lats_south),plotarray_south,color='red'  ,linestyle='dashdot')

        axs[0,0].tick_params(axis='both',labelsize=tick_font_size)    
        axs[0,0].set_ylabel(r'\textbf{Temperature (K)}',size = axis_text_size)
        axs[0,0].set_xlabel(r'\textbf{Latitude ($^{/circ}$)}',size = axis_text_size)
        axs[0,0].set_xlim(right=90.0)
        axs[0,0].set_xlim(left= 0.0)
        
        plotarray = NT[:nLats,iAlt]

        plotarray_north = np.ma.masked_where(lats <= 0.0, plotarray)
        lats_north      = np.ma.masked_where(lats <= 0.0, lats)

        plotarray_south = np.ma.masked_where(lats >= 0.0, plotarray)
        lats_south      = np.ma.masked_where(lats >= 0.0, lats)

     #   axs[0,1].plot(lats,plotarray)
        axs[0,1].plot(    lats_north ,plotarray_north,color='blue' ,linestyle='dashed')
        axs[0,1].plot(abs(lats_south),plotarray_south,color='red'  ,linestyle='dashdot')
        axs[0,1].tick_params(axis='both',labelsize=tick_font_size)    
        axs[0,1].set_ylabel(r'\textbf{NDensity }',size = axis_text_size)
        axs[0,1].set_xlabel(r'\textbf{Latitude ($^{/circ}$)}',size = axis_text_size)
       # axs[0,1].set_xlim(right=90.0)
       # axs[0,1].set_xlim(left=-90.0)
        axs[0,1].set_xlim(right=90.0)
        axs[0,1].set_xlim(left= 0.0)
        
        
        plotarray = U[:nLats,iAlt]
        plotarray_north = np.ma.masked_where(lats <= 0.0, plotarray)
        lats_north      = np.ma.masked_where(lats <= 0.0, lats)

        plotarray_south = np.ma.masked_where(lats >= 0.0, plotarray)
        lats_south      = np.ma.masked_where(lats >= 0.0, lats)

        #axs[1,0].plot(lats,plotarray)
        axs[1,0].plot(    lats_north ,plotarray_north,color='blue' ,linestyle='dashed')
        axs[1,0].plot(abs(lats_south),plotarray_south,color='red'  ,linestyle='dashdot')
     
        axs[1,0].tick_params(axis='both',labelsize=tick_font_size)    
        axs[1,0].set_ylabel(r'\textbf{Velocity East (m/s)}',size = axis_text_size)
        axs[1,0].set_xlabel(r'\textbf{Latitude ($^{/circ}$)}',size = axis_text_size)
        #axs[1,0].set_xlim(right=90.0)
        #axs[1,0].set_xlim(left=-90.0)
        axs[1,0].set_xlim(right=90.0)
        axs[1,0].set_xlim(left= 0.0)

        plotarray = V[:nLats,iAlt]
        plotarray_north = np.ma.masked_where(lats <= 0.0, plotarray)
        lats_north      = np.ma.masked_where(lats <= 0.0, lats)

        plotarray_south = np.ma.masked_where(lats >= 0.0, plotarray)
        lats_south      = np.ma.masked_where(lats >= 0.0, lats)
        
       # axs[1,1].plot(lats,plotarray)
        axs[1,1].plot(    lats_north ,plotarray_north,color='blue' ,linestyle='dashed')
        axs[1,1].plot(abs(lats_south),plotarray_south,color='red'  ,linestyle='dashdot')
     
        axs[1,1].tick_params(axis='both',labelsize=tick_font_size)    
        axs[1,1].set_ylabel(r'\textbf{Velocity North (m/s)}',size = axis_text_size)
        axs[1,1].set_xlabel(r'\textbf{Latitude ($^{/circ}$)}',size = axis_text_size)
        #axs[1,1].set_xlim(right=90.0)
        #axs[1,1].set_xlim(left=-90.0)
        axs[1,1].set_xlim(right=90.0)
        axs[1,1].set_xlim(left= 0.0)

        #plt.subplots_adjust(hspace=0.1, right = 1.030)
        #plt.subplots_adjust(bottom = 0.12, right = 1.0, top = 0.90, left=0.12)
        pdf.savefig()
        plt.close()
        AltSliceIndex = AltSliceIndex + 1

print('Finished Output to File: '+ pdfile)
LonSliceIndex = 0





#
#
#

