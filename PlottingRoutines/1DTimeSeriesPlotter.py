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

## for debugging
import pdb


# Define a short routine for asking inputs
def ask_user(question):
    check = str(input(question+' (y/n): ')).lower().strip()
    if check[:1] == 'y':
       return True
    if check[:1] == 'n':
       return False
    else:
       print('Invalid Input')
       return ask_user("LMAO.. No really, Plot Log?")


# First, we find all 1DALL binary files in the directory and then
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
#print('===> Searching Working Directory for GITM 3DALL Files. \n')
filelist = []  # initialize a filelist array
for file in glob.glob("1DALL*.bin"):
    filelist.append(file)
    
# Sort the Files by date
filelist.sort()
if len(filelist) > 0:
   print('===> Identified ', len(filelist), ' Gitm Binary Files.\n')
   testfile = []
   for file in filelist:
       testfile.append(file)
       #print(file)
else:
    print("!!! Couldn't find any 3DALL*.bin GITM files !!!!")
    raise Exception()
    
gdata = gitm.GitmBin(testfile[len(filelist) - 1])
OrignVars = gdata.attrs['nVars']
# Pull out the keys (Variable Names) to use
Variables = gdata.keys()
OrigVariableNames = list(Variables)

VariableNames = OrigVariableNames[4:OrignVars]
#print(VariableNames)

# nTimes is the number of files that we have
# Determines how many timestamps we average over
nTimes = len(filelist)

# Determine how many Altitude points are present
OrignAlts = gdata.attrs['nAlt']
alts = gdata['Altitude' ][0     ,0     ,0:OrignAlts]/1000.0    # Alts in km
nAlts = len(alts)


# Next, we create the Time-Averaged Mean of our Simulation
# Set the shape of our array
# Note that we ignore the first 4 variables
#
nVars = OrignVars - 4
array_shape = (nTimes, nVars, nAlts)
VariableProfile = np.zeros(array_shape,order='F')

for iFile in range(nTimes):
    gdata = gitm.GitmBin(testfile[iFile])
    for iVar in range(nVars):
        KeySelect = VariableNames[iVar]
        VariableProfile[iFile,iVar,0:nAlts] = \
              gdata[KeySelect][0,0,0:nAlts]

# Now VariableTimeSeries has our data as [nTimes, nVars, nAlts]
# Next, we average over time

# Mean Variable Profiles is now the time-averged profile
# size is (nVars, nAlts)

# Our Logical Flag for Plotting
KeepPlotting = True
PI = np.pi

# Some Universal Plot Settings
axis_text_size = 18
title_text_size = 28
tick_font_size = 13

# GET USER INPUT

#pdb.set_trace()
while KeepPlotting == True :

   for iVar in range(nVars):
       print('[',iVar,']   ',VariableNames[iVar])
   UserSelectionString = input("Pick a Variable to Plot [Enter Number]:  ")
   iUser_ = int(UserSelectionString) 
   UserKey = VariableNames[iUser_]
   print(' \n')
   print('You want to Plot: ',UserKey,'\n')

   # Now we have the Time Series of the Selected Variable
   VariableString = UserKey

   question = " Plot Logarithmically? "
   plotlog = ask_user(question)
 
   TimeSeriesVariable = VariableProfile[0:nTimes,iUser_,0:nAlts]

   pdfile = '1D_TimeSeries_' + UserKey.strip() + '.pdf'
   TitleString = 'Time Series Mean ' + UserKey.strip() 
   if plotlog == True:

       X = np.linspace(0,1,len(lats))
       Y = np.linspace(0,1,len(alts))
       for i in range(0,len(lats)):
           X[i] = lats[i]
       for j in range(0,len(alts)):
           Y[j] = alts[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
       array_shape = (len(lats), len(alts))
       X2D = np.zeros(array_shape,order='F')
       Y2D = np.zeros(array_shape,order='F')
       Z2D = np.zeros(array_shape,order='F')
       for i in range(0,nLats):
           for j in range(0,nAlts):
               X2D[i,j] = X[i]
               Y2D[i,j] = Y[j]
        
       for i in range(0,len(lats)):
           for j in range (0, len(alts)):
               Z2D[i,j]= Temperature[i,j]

#      fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,8.5),sharex=True)
#      line1, = ax.plot(TimeAvgVariable, alts)
#      ax.tick_params(axis='both',labelsize=tick_font_size)    
#      ax.set_title(TitleString,size = title_text_size)   
#      ax.set_xlabel(''+UserKey,size = axis_text_size)
#      ax.set_ylabel('Altitude (km)',size = axis_text_size)
#      plt.xscale("log")

#      plt.show()
#      question = " Save to File? "
#      savefile = ask_user(question)
#      if savefile == True:
#         plt.savefig(pdfile)
#      plt.close()



   else:
       
       X = np.linspace(0,1,nTimes)
       Y = np.linspace(0,1,len(alts))
       for i in range(0,nTimes):
           X[i] = i
       for j in range(0,len(alts)):
           Y[j] = alts[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
       array_shape = (nTimes, len(alts))
       X2D = np.zeros(array_shape,order='F')
       Y2D = np.zeros(array_shape,order='F')
       Z2D = np.zeros(array_shape,order='F')
       for i in range(0,nTimes):
           for j in range(0,nAlts):
               X2D[i,j] = X[i]
               Y2D[i,j] = Y[j]
        
       for i in range(0,nTimes):
           for j in range (0, len(alts)):
               Z2D[i,j]= TimeSeriesVariable[i,j]
       fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,5.0),sharex=True)
       bw1 = ax.contour(X2D,Y2D,Z2D,10,colors='k',linewidths=1.5)
       cpf = ax.contourf(X2D,Y2D,Z2D,20,cmap=cmocean.cm.thermal)
       # Set axis text size and the labels
       TitleString = '' + UserKey + ' Time Series'
       ax.annotate(TitleString,xy=(0.1,0.95),xycoords='figure fraction',fontsize=title_text_size)
       ax.tick_params(axis='both',labelsize=tick_font_size)
       ax.set_ylabel('Altitude (km)',size = axis_text_size)
        
       # Tell the colorbar where to exist
       # fig.colorbar(subplotname, axis name)
       cbar = fig.colorbar(cpf,ax=ax)   
       cbar.ax.tick_params(labelsize='12')
        
       ax.set_xlabel('Simulation Time (Earth Days)',size = axis_text_size)

       plt.show()
       question = " Save to File? "
       savefile = ask_user(question)
       if savefile == True:
          plt.savefig(pdfile)
       plt.close()
        #fig.savefig('GITM_Sample_Contour_Plasma.eps')
       
 #     fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,8.5),sharex=True)
 #     line1, = ax.plot(TimeAvgVariable, alts)
 #     ax.tick_params(axis='both',labelsize=tick_font_size)    
 #     ax.set_title(TitleString,size = title_text_size)   
 #     ax.set_xlabel(''+UserKey,size = axis_text_size)
 #     ax.set_ylabel('Altitude (km)',size = axis_text_size)
#
#      plt.show()
#      question = " Save to File? "
 #     savefile = ask_user(question)
#      if savefile == True:
#         plt.savefig(pdfile)
#     plt.close()

   question = " Plot Another Variable? "
   KeepPlotting = ask_user(question)

print('GITM_Time_Series 1D PLOTTER:  Quitting Plot Routine')
