#!/usr/bin/env python
"""
Created on Mon Oct 16 14:35:06 2017

This python script tranform the Particles-###.vtk files in to csv files

extractInfoVTK(vtkfilename) return the velocities, the volumes and points of each MPM point
readconfigFile(configFile) return a disctionary with 18 firs lines of the config file 

@author: Gauthier
"""
#%%
import numpy as np  
np.array([10,10])
import json
import sys, os
import subprocess
driveFolder='/scratch/garousse/'
# driveFolder='/media/gauthier/Gauthier_Backup/'
sys.path.append(driveFolder+'TAF/TAF_EPFL/current_work/OPyF-Project/github/opyFlow/')
sys.path.append(driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/python-essentials/Essentials OpenCV')
#import trackandinterpolate
 
 #the path where all the videos are

vidPath=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Video_src' 
 # the genarated dictionnary-json path
JSONpath=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Video_src/dictExp.json'

    # the d6 soft path
d6Path=driveFolder+'TAF/TAF_inria/Sand6/epfl_lhe_2d_and_3d/build_fast'
d6Path='/media/gauthier/Data-Gauthier/programs/gitLab/sand6/build'
d6Path='/scratch/garousse/TAF/TAF_inria/INRIA_current_work/GitLab/sand6/build'
#d6Path=driveFolder+'TAF/TAF_inria/GitLab/sand6cohesive/build_julien'
#d6Path='/home/gauthier/programs/epfl_lhe/build2d'
d6OutFolder='out'
out_Opyf=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs_opyf/'


os.chdir(d6Path)
    #add d6py module pythonpath
sys.path.append(d6Path+'/../python')

import d6py
#%%

runName='test_sans_hyst_I0_0.005_deltamuS_0.07_sans_lit_en_aval_condtion_au_bord_bloc_low_mu_Edwards_exp'

d6OutFolder=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs/Tests/'+runName
d6py.mkdir2(d6OutFolder) 

delta_mu=0.
I0=0.3
I0_start=0.005
# delta_mu_start=np.tan(25.3*np.pi/180)-np.tan(21.3*np.pi/180)
delta_mu_start=0
muRigid=0.
P0=1.
    
Lmod=0.6
Hmod=0.1
import time
t=time.time()
slope=np.arctan(0.45)*180/np.pi

configFilein=d6Path+'/../scenes/hyteresis_test.3d.LHE.Door.conf'
newConfigFile=d6OutFolder+'/hyteresis_test_'+runName+'m.conf'

d6py.modifyConfigFile(configFilein,newConfigFile,'box',[Lmod, 0.06,Hmod])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'fps',[60])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'gravity',[+9.81*np.sin(slope*np.pi/180), 0.000,-9.81*np.cos(slope*np.pi/180)])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'nFrames',[200])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'randomize',[0])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'substeps',[1])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'grainDiameter',[0.0005])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'mu',[np.tan(21.3*np.pi/180)])
# d6py.modifyConfigFile(newConfigFile,newConfigFile,'mu',[np.tan(25.3*np.pi/180)])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'muRigid',[muRigid]) 
d6py.modifyConfigFile(newConfigFile,newConfigFile,'delta_mu',[delta_mu])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'delta_mu_start',[delta_mu_start])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'I0_start',[I0_start])
d6py.modifyConfigFile(newConfigFile,newConfigFile,'P0',[P0])

d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','rb_plane_test')

TypicalLength=0.01
d6py.modifyConfigFile(newConfigFile,newConfigFile,'res',[Lmod//TypicalLength,0.06//TypicalLength,np.round(Hmod/TypicalLength/10)*10]) #pour avoir un réolution divisible par 10 selon Y
d6py.modifyConfigFile(newConfigFile,newConfigFile,'I0',[I0]) 
  
#load the final config file dictionnary    
    
dConfigmod=d6py.readConfigFile(newConfigFile)

d6py.d6run(d6OutFolder,newConfigFile)

d6py.d62vtk(d6OutFolder,allF=True,particles=True)


    
    
    
    

