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


#%
# Load the contents from the file, which creates a new dictionary with all experimental infos
in_file = open(JSONpath,"r")
dictExp = json.load(in_file)   


#generate the list of avi files
lExp=[]
for d in dictExp:
    lExp.append(d)        
lExp=np.sort(lExp) 

import time
t=time.time()

door='with'

for j in range(7,8 ): 
#    plt.close('all')
#for j in range(0,8):
    sE=lExp[j] #Selected exeperiment
    sdictE=dictExp[sE]
    fracH=0.8
    
    if sdictE['camType']=='Phantom':
        sdictE['L']=sdictE['L']+0.1
        sdictE['Ltot']=sdictE['Ltot']+0.1
    #We Consider a 1 meter simu
    L=sdictE['L']    
    Lmod=sdictE['Ltot']
    nFrames=sdictE['nFrames']
    nFrames=1
    
    if (sdictE['camType']=='BW') & (sdictE['Slope']==15. or sdictE['Slope']==20.):
        Lmod=sdictE['Ltot']+0.5

    delta_mu=0.2
    I0=0.3
    I0_start=0.003
    delta_mu_start=0.07
    muRigid=0.18
    P0=1.
    
    Hmod=(sdictE['H']+0.005)/fracH
    
    if door=='with':
        ts=0
    else:
        ts=20
    substeps=1
    sdictE['delta_mu']=delta_mu

    mu=sdictE['mu']
    prop=''
    if prop=='low':
        mu=sdictE['mu']-0.05
    
    if prop=='inscrit':
        mu=mu/(1+1/3*mu**2)

    prop='Test_P0'   
    
    if door=='with':
        runName=str('Run_'+format(j,'02.0f')+'_3D_Door_muRigid='+str(muRigid)+'_'+sdictE['grainType']+'_Slope='+format(sdictE['Slope'],'.0f')+'deg_H='+format(sdictE['H'],'.3f')+'m_L='+format(sdictE['L'],'.3f')+'_delta_mu='+format(delta_mu,'.3f')+'_substeps_'+str(substeps)+'_fracH='+str(fracH)+'_I0_start='+format(I0_start,'.4f')+'_delta_mu_start='+format(delta_mu_start,'.4f')+'_P0='+format(P0,'.4f')+prop)    
    else:
        runName=str('Run_'+format(j,'02.0f')+'_3D_no_Door_start_at_0.13_s_'+sdictE['grainType']+'_Slope='+format(sdictE['Slope'],'.0f')+'deg_H='+format(sdictE['H'],'.3f')+'m_L='+format(sdictE['L'],'.3f')+'_delta_mu='+format(delta_mu,'.3f')+'_substeps_'+str(substeps)+'_fracH='+str(fracH)+'_I0_start='+format(I0_start,'.4f')+'_delta_mu_start='+format(delta_mu_start,'.4f')+'_P0='+format(P0,'.4f')+prop)      
       
     
    d6OutFolder=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs/Tests/'+runName
    d6py.mkdir2(d6OutFolder) 
    
    newConfigFile=d6OutFolder+'/collapse.3d_'+runName+'m.conf'
    configFilein=d6Path+'/../scenes/collapse.3d.LHE.Door.conf'

    d6py.modifyConfigFile(configFilein,newConfigFile,'box',[Lmod, 0.06,Hmod])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'fps',[150])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'gravity',[+9.81*np.sin(sdictE['Slope']*np.pi/180), 0,-9.81*np.cos(sdictE['Slope']*np.pi/180)])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'nFrames',[nFrames+ts])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'randomize',[0])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'substeps',[substeps])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'grainDiameter',[sdictE['grainDiameter']])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'mu',[mu])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'muRigid',[muRigid]) #mu Door
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'delta_mu',[delta_mu])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'delta_mu_start',[delta_mu_start])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'I0_start',[I0_start])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'P0',[P0])
    if door=='with':
        if j<=3:
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.06 veldoor:0.8 ts:'+format(ts,'1.0f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L))
        else:
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.13 veldoor:0.7 ts:'+format(ts,'1.0f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L))
    else:
        d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.0001 veldoor:100 ts:'+format(ts,'1.0f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L))

    
    TypicalLength=0.005
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'res',[Lmod//TypicalLength,0.06//TypicalLength,np.round(Hmod/TypicalLength/10)*10]) #pour avoir un réolution divisible par 10 selon Y
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'I0',[I0]) 
      
    #load the final config file dictionnary    
        
    dConfigmod=d6py.readConfigFile(newConfigFile)
    #% 

    ################# Execute Different sytem commands ###############     
     
    #proc =subprocess.Popen('rm -rf  ' + d6OutFolder +'/*' , shell=True,stdout=subprocess.PIPE)
    #(out, err) = proc.communicate()    
    # d6py.d6run(d6OutFolder,newConfigFile)
    
    d6py.d62vtk(d6OutFolder,allF=True,particles=True)

    # cmd = str('./apps/d6 ' +d6OutFolder+ ' -i '+ newConfigFile)
    
    # proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # (out, err) = proc.communicate()
    # print("program output:", out)
    
    
    # if out=='':
    #     print('No sand6 simulation were executed - run aborted')
    #     sys.exit()  
    # else:
    #     print("program output:", out) 
    
    
    
    # cmd = str('./apps/d62vtk '+d6OutFolder+ ' -a -p')
    # proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    # (out, err) = proc.communicate()
     
    
    # if out=='':
    #     print('No vtk files generated - run aborted')
    #     sys.exit() 
    # else:
    #     print("program output:", out)

    
    
    
    

