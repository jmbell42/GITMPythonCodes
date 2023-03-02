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
OrignVars = gdata.attrs['nVars']
Variables = gdata.keys()
OrigVariableNames = list(Variables)
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


# These are the indices of the Altitudes that we want to pull out.

AltString = list()


pdfile = 'ConstantAltitude_DensityContours.pdf'
with PdfPages(pdfile) as pdf:
    #for i in AltitudeSlices:
    #for i in range(nAlts):
    for i in range(5):

        altpick = alts[i]
        newstring = str(altpick)
        splitstring = newstring.split('.')
        AltString.append(splitstring[0])

    AltSliceIndex = 0
    for iAlt in range(nAlts):
        print('Accessing iAlt %i',iAlt)
        iAlt = min(iAlt,nAlts-1)
#       Temperature = gdata['Temperature'   ][:nLons,:nLats,iAlt]   # Pull out Temperature
#       U           = gdata['V!Dn!N (east)' ][:nLons,:nLats,iAlt]
#       V           = gdata['V!Dn!N (north)'][:nLons,:nLats,iAlt]
#       W           = gdata['V!Dn!N (up)'][:nLons,:nLats,iAlt]
#        Rho         = gdata['Rho'][:nLons,:nLats,iAlt]
        nN2 = gdata['N2'   ][:nLons,:nLats,iAlt]   # Pull out Temperature
        nCH4 = gdata['CH4' ][:nLons,:nLats,iAlt]   # Pull out Temperature
        nH2 = gdata['H2'   ][:nLons,:nLats,iAlt]   # Pull out Temperature
        nNT = nN2 + nCH4 + nH2

        # Create local numpy arrays for plotting
        X = np.linspace(0,1,len(lons))
        Y = np.linspace(0,1,len(lats))
        for i in range(0,len(lons)):
            X[i] = lons[i]
        for j in range(0,len(lats)):
            Y[j] = lats[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
        array_shape = (len(lons), len(lats))
        X2D = np.zeros(array_shape,order='F')
        Y2D = np.zeros(array_shape,order='F')
        Z2D = np.zeros(array_shape,order='F')
        for i in range(0,nLons):
            for j in range(0,nLats):
                X2D[i,j] = X[i]
                Y2D[i,j] = Y[j]
        for i in range(0,len(lons)):
            for j in range (0, len(lats)):
                Z2D[i,j]= nNT[i,j]
                #Z2D[i,j]= Temperature[i,j]
        # Some plot parameters
        axis_text_size = 18
        title_text_size = 28
        tick_font_size = 13
        projection = ccrs.PlateCarree()
        # Figsize = (width, height) in inches
        # subplots (nrows, ncolumns, figure_settings)
        # Instantiate a figure container.  Fig is the whole figure
        # ax is each sub figure.
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,11),sharex=True)
        bw1 = ax.contour(X2D,Y2D,Z2D,10,colors='k',linewidths=1.5)
        cpf = ax.contourf(X2D,Y2D,Z2D,20,cmap=cmocean.cm.balance)
        # Set axis text size and the labels
        TitleString = 'Altitude = ' + AltString[AltSliceIndex] + ' km'
        ax.tick_params(axis='both',labelsize=tick_font_size)
        ax.set_title(TitleString,size = title_text_size)
        ax.set_xlabel(r'\textbf{Longitude (deg)}',size = axis_text_size)
        ax.set_ylabel(r'\textbf{Latitude (deg)}',size = axis_text_size)
        ax.set_xlim(right=360.0)
        ax.set_xlim(left=0.0)
        ax.set_ylim(top=  90.0)
        ax.set_ylim(bottom=-90.0)
        
        # Tell the colorbar where to exist
        # fig.colorbar(subplotname, axis name)
        cbar = fig.colorbar(cpf,ax=ax,pad=0.02,format='%.4e')   
        cbar.ax.tick_params(labelsize='14')
        cbar.set_label(r'Density ',size='large',weight='bold')
    
        plt.subplots_adjust(hspace=0.12, right = 0.95, top = 0.90, left = 0.12)
        pdf.savefig()
        plt.close()

        AltSliceIndex = AltSliceIndex + 1
print('Finished Output to File: '+ pdfile)




