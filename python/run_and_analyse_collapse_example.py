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

configFileOrigin='../scenes/collapse.3d.LHE.Door.conf'
d6OutFolder='../build/out' # create the folder where you want 

configFileModified='../scenes/collapse.3d.LHE_m.Door.conf'
shutil.copy(configFileOrigin, configFileModified) 

slope=10
d6py.modifyMultipleConf(configFileModified,viscosity=[0.0], volMass=[1500], 
                    phiMax=[1], fps=[15], 
                    gravity=[+9.81*np.sin(slope*np.pi/180), 0,-9.81*np.cos(slope*np.pi/180)],
                    nFrames=[2],randomize=[1],
                    substeps=[10], grainDiameter=[0.0005],
                    mu=[0], muRigid=[0.23], delta_mu=[0], delta_mu_start=[0],
                    I0_start=[0],useInfNorm=[1], nSamples=[2],)


pathlib.Path(d6OutFolder).mkdir(parents=True, exist_ok=True)

d6run(d6OutFolder,configFileModified)


#%%



