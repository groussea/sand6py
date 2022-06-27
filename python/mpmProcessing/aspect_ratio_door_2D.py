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
import json
import sys, os

fileDir = os.path.dirname(os.path.abspath(__file__))

print(fileDir)

d6Path=fileDir+'/../../build-2d-p39'

mainOutFolder=d6Path+'/out'

JSONpath=fileDir+'/../Granular_Collapses_Experimental_Informations.json'


os.chdir(d6Path)
    #add d6py module pythonpath
sys.path.append(d6Path+'/../python')
import d6py
from d6py.d6python2D import * # python must be reload if d6python2D was imported

mainOutFolder='/media/gauthier/DataSSD/sand6_out/2d/'

#%%
def rund6py(sdictE,**args):
    fracH=0.8

    H0 = args.get('H0', 0.12)
    Hmod = H0/fracH  # AR1_ experimental
    # Hmod=0.12/fracH #AR1_ experimental
    L0 = args.get('L0', H0)
    L = L0   # to add a wall on the left
    
    Lmod = np.round((H0/L0)*L0+3*L0,2) 
    nFrames = 35
    
        
    if (sdictE['camType']=='BW') & (sdictE['Slope']==15. or sdictE['Slope']==20.):
        Lmod=sdictE['Ltot']+0.8

    delta_mu=args.get('delta_mu',0)
    I0=args.get('I0',0.3)
    I0_start=args.get('I0_start',0.000)
    delta_mu_start=args.get('delta_mu_start',0.00)
    muRigid=args.get('muRigid',0.18)
    mu=args.get('mu',sdictE['mu'])
    
    P0=1.
    
    door=args.get('door','with')
    if door=='with':
        ts=0.2
    else:
        ts = 0.3 # Attention à modifier car en secondes
        

    substeps=args.get('substeps',40)
    resZ=args.get('resZ',60)

    prop=args.get('prop','test')
    rand = args.get('rand', 0)
    
    if door=='with':
        runName=str('Run_'+format(j,'02.0f')+'_2D_Door_mu='+str(mu)+'_H_'+format(Hmod*100,'.0f')+'_L_'+format(L*100, '.0f')+'_cm_Slope='+format(sdictE['Slope'],'.0f')+'_fracH='+'_R_'+format(H0/L0, '.1f')+'_'+prop)    
    else:
        runName=str('Run_'+format(j,'02.0f')+'_2D_no_Door_mu='+str(mu)+'_H_'+format(Hmod*100,'.2f')+'_L_'+format(L, '.2f')+'cm_Slope='+format(sdictE['Slope'],'.0f')+'_fracH='+'_R_'+format(H0/L0, '.1f')+'_'+prop)   
       
     
    d6OutFolder=mainOutFolder+runName
    d6py.mkdir2(d6OutFolder) 
    
    newConfigFile=d6OutFolder+'/collapse.2d_'+runName+'m.conf'
    configFilein=d6Path+'/../scenes/collapse.2d.LHE.Door.conf'

    d6py.modifyConfigFile(configFilein,newConfigFile,'box',[Lmod, Hmod])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'fps',[15])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'volMass',[sdictE['density']])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'gravity',[+9.81*np.sin(sdictE['Slope']*np.pi/180), -9.81*np.cos(sdictE['Slope']*np.pi/180)])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'nFrames',[nFrames+ts*args.get('fps',15)])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'randomize',[rand])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'substeps',[substeps])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'grainDiameter',[sdictE['grainDiameter']])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'mu',[mu])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'muRigid',[muRigid]) #mu Door
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'delta_mu',[delta_mu])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'delta_mu_start',[delta_mu_start])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'I0_start',[I0_start])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'P0',[P0])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'nSamples',[args.get('nSamples',3)])
    if door=='with':
        if j<=3:
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.06 veldoor:0.8 ts:'+format(ts,'1.2f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L) +' hbed:'+str(0.0) )
        else:
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.13 veldoor:0.7 ts:'+format(ts,'1.2f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L)+' hbed:'+str(0.0) )
    else:
        d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.0001 veldoor:100 ts:'+format(ts,'1.2f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L)+' hbed:'+str(0.0) )

    
    TypicalLength=0.01
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'res',[int(Lmod/TypicalLength),int(Hmod/TypicalLength)]) #pour avoir un réolution divisible par 10 selon Y
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'I0',[I0]) 
      
    #load the final config file dictionnary    
        
    # dConfigmod=d6py.readConfigFile(newConfigFile)
    #% 

    d6run(d6OutFolder,newConfigFile)
    
    # d6py.d62vtk(d6OutFolder,allF=True,particles=True)



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

# for j in range(0,1 ):  
#     sE=lExp[j] #Selected exeperiment
#     sdictE=dictExp[sE]
#     for dmu in [0.05,-0.05]:
#         for s,p in zip([40],['fin']):
#             rund6py(sdictE, delta_mu=0.,rand=0,mu=np.round(sdictE['mu']+dmu,2),substeps=s,prop=p,muRigid=0.,nSamples=3)
            
for j in [7]:  
    sE=lExp[j] #Selected exeperiment
    sdictE=dictExp[sE]
    for m in [0.44, 0.45, 0.46, 0.48, 0.5,0.54, 0.58]:
        # rund6py(sdictE,delta_mu=0.,rand=0,mu=m,substeps=120,prop='aspect_ratio_1',muRigid=0.,nSamples=15,door='with')
        L0=0.12
        for R in [1]:
            rund6py(sdictE, delta_mu=0., mu=m, prop='hbed_0', fps=15, nFrames=25, nSamples=15, I0_start=0.005, delta_mu_start=0.0, visc=0.0, rand=0, substeps=120, I0=0.279, delta_x=0.01, H0=R*L0, L0=L0)
        # rund6py(sdictE,delta_mu=0.,rand=0,mu=0.65,substeps=s,prop=p,muRigid=0.,nSamples=15,door='with')


# for j in range(7,8 ):  
#     sE=lExp[j] #Selected exeperiment
#     sdictE=dictExp[sE]
#     for dmu in [0.05,-0.05]:
#         for s,p in zip([20],['rand-longer-higher-res']):
#             rund6py(sdictE,delta_mu=0.,rand=0,mu=np.round(sdictE['mu']+dmu,2),substeps=s,prop=p,muRigid=0.,nSamples=3)
# for j in range(7,9 ):  
#     sE=lExp[j] #Selected exeperiment
#     sdictE=dictExp[sE]
#     for dmu in [0.05,-0.05]:
#         for s,p in zip([20],['rand-longer-higher-res']):
#             rund6py(sdictE,delta_mu=0.,rand=1,mu=np.round(sdictE['mu']+dmu,2),substeps=s,prop=p,muRigid=0.)

 #%%   
    #    
    
#%%




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

    
    
    
    

