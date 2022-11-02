#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 2  2022

@author: Gauthier Rousseau
"""

#%%
import d6py
import pathlib
import shutil
import numpy as np
from d6py.d6python3D import * # python must be reload if d6python2D was imported
# from d6py.d6python2D import * # if 2D

configFileOrigin='../scenes/collapse.3d.LHE.Door.conf'
d6OutFolder='../build/out' # create the folder where you want 

configFileModified='../scenes/collapse.3d.LHE_m.Door.conf'
shutil.copy(configFileOrigin, configFileModified) 

slope=0
d6py.modifyMultipleConf(configFileModified,res=[100, 8,	20],
                        box=[1.0, 0.1/0.8, 0.10	],
                        viscosity=[0.0], volMass=[1500], 
                        phiMax=[1], fps=[15], 
                        gravity=[+9.81*np.sin(slope*np.pi/180), 0,-9.81*np.cos(slope*np.pi/180)],
                        nFrames=[20],randomize=[0],
                        substeps=[10], grainDiameter=[0.0005],
                        mu=[0.44], muRigid=[0.23], delta_mu=[0], 
                        delta_mu_start=[0],
                        I0_start=[0],useInfNorm=[0], nSamples=[2],)

L=0.1 #length of the column
ts=0.2 #time to wait before opening gate
d6py.modifyConfigFile(configFileModified,configFileModified,'scenario','collapselhedoor taudoor:0.0001 veldoor:100 ts:'+format(ts,'1.2f')+' frac_h:'+format(0.8,'1.1f')+' column_length:'+str(L))

pathlib.Path(d6OutFolder).mkdir(parents=True, exist_ok=True)

d6run(d6OutFolder,configFileModified)
#running may take while
#run the following code with another python session for analysing during the process

#%%

import matplotlib.pyplot as plt

run=d6py.NumericalRun(d6OutFolder)
path=pathlib.Path(d6OutFolder+'/vtk') #remove vtk folder if it contains previous simulation
shutil.rmtree(path)
run.loadVTK(4) #generate the vtk file if it does not exist

fig, ax = plt.subplots(1,1)

# run.plotVelocity(ax)

# run.plotContour(ax,levels=[0.5])
run.opyfPointCloudColoredScatter(ax)
#plot the velocity


#%% more elaborated analysis




