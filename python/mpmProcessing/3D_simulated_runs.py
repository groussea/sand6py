
#!/usr/bin/env python
"""
Created on Mon Oct 16 14:35:06 2017


Python script to run the several cases of JFM
readconfigFile(configFile) return a dictionary 
@author: Gauthier
"""
#%%

from re import A
import numpy as np  
import json
import sys, os
sys.path.append('./..')
import d6py
from d6py.Tools import *
fileDir = os.path.dirname(os.path.abspath(__file__))

print(fileDir)

d6Path=fileDir+'/../../sand6/build3d'

mainOutFolder=d6Path+'/out'

JSONpath=fileDir+'/../Granular_Collapses_Experimental_Informations.json'


from d6py.d6python3D import * # python must be reload if d6python2D was imported

#% by default, mainOutFolder is in the build folder but it is highly recommended to set your own mainOutPut folder since outputs are generally large
# mainOutFolder=d6Path+'/out'
# mainOutFolder='/media/gauthier/DataSSD/sand6_out/'
mainOutFolder=fileDir+'/../../collapses_simulations/constant_mu_avalanche'
# mainOutFolder = '/media/gauthier/SSD500/sand6_out/'
# mainOutFolder='/home/gauthier/sorties_sand6/'
d6py.mkdir2(mainOutFolder) 
#%%

def modifyMultipleConf(configFile,**args):
    for a in args:
        d6py.modifyConfigFile(configFile,configFile,a,args[a])
    

fignames=['G00', 'G05', 'G10', 'G15', 'B00', 'B05', 'B10', 'B15', 'B20']
#%%
def rund6py(sdictE,**args):
    fracH=0.8
    if sdictE['camType']=='Phantom':
        sdictE['L']=sdictE['L']+0.1
        sdictE['Ltot']=sdictE['Ltot']+0.1

        
    L=sdictE['L']+0.01   # to add a wall on the left 
    Lmod=sdictE['Ltot']
    nFrames=args.get('nFrames',int(sdictE['nFrames']/10))
        
    if (sdictE['camType']=='BW') & (sdictE['Slope']==15.):
        Lmod = sdictE['Ltot'] + 0.6
        
    if (sdictE['camType']=='BW') & (sdictE['Slope']==20.):
        Lmod = sdictE['Ltot'] + 1.2
        
    delta_mu=args.get('delta_mu',0)
    I0=args.get('I0',0.279)
    I0_start=args.get('I0_start',0.000)
    delta_mu_start=args.get('delta_mu_start',0.00)
    muRigid=args.get('muRigid',0.18)
    mu=args.get('mu',sdictE['mu'])
    
    P0=1.
    
    Hmod=(np.round(sdictE['H'],3))/fracH
    door=args.get('door','with')
    if door=='with':
        ts=0.2
    else:
        ts = 2 # to change with the fps
        

    substeps=args.get('substeps',20)
    resZ=args.get('resZ',30)
    wsw=args.get('wsw',0.01)
    muDoor = args.get('mudoor',-99)
    prop=args.get('prop','test')
    
    visc=args.get('visc',0.0)
    
    rand = args.get('rand', 1)
    W = args.get('W', 0.08)


    runName=fignames[j]+'_3D_mu='+str(mu)+'_muRigid='+str(muRigid)+'_W_'+format(W*100,'.1f')+'cm_'+sdictE['grainType']+'_Slope='+format(sdictE['Slope'],'.0f')+'deg_delta_mu='+format(delta_mu,'.2f')+'_substeps_'+str(substeps)+'_I0_start='+format(I0_start,'.4f')+'_delta_mu_start='+format(delta_mu_start,'.2f')+prop    


    d6OutFolder=mainOutFolder+'/'+runName
    d6py.mkdir2(d6OutFolder) 
    
    newConfigFile=d6OutFolder+'/collapse.3d_'+runName+'.conf'
    configFilein=d6Path+'/../scenes/collapse.3d.LHE.Door.conf'
    
    
    
    d6py.modifyConfigFile(configFilein,newConfigFile,'box',[Lmod, W,Hmod])
    
    modifyMultipleConf(newConfigFile,viscosity=[visc], volMass=[sdictE['density']], 
                       phiMax=[1], fps=[args.get('fps',15)], 
                       gravity=[+9.81*np.sin(sdictE['Slope']*np.pi/180), 0,-9.81*np.cos(sdictE['Slope']*np.pi/180)],
                       nFrames=[nFrames+ts*args.get('fps',15)],randomize=[rand],
                       substeps=[substeps], grainDiameter=[sdictE['grainDiameter']],
                       mu=[mu], muRigid=[muRigid], delta_mu=[delta_mu], delta_mu_start=[delta_mu_start],
                       I0_start=[I0_start],useInfNorm=[1], P0=[P0], nSamples=[args.get('nSamples',3)],)
    

    if door=='with':
        if j<=3:
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.06 veldoor:0.8 ts:'+format(ts,'1.2f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L)+' wsw:'+str(wsw) +' mud:'+str(muDoor))
        else:
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.13 veldoor:0.7 ts:'+format(ts,'1.2f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L)+' wsw:'+str(wsw)+' mud:'+str(muDoor))
    else:
        d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.0001 veldoor:100 ts:'+format(ts,'1.2f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L))

    d6py.modifyConfigFile(newConfigFile,newConfigFile,'boundary','top:slip left:slip right:stick front:slip back:slip bottom:stick')
    
    delta_x=args.get('delta_x',0.01)
    delta_x=L/np.round(L/0.01)
    delta_y=args.get('delta_y',0.01)
    
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'res',[int(Lmod/delta_x),int(W/delta_y),resZ])
    #  pour avoir un réolution divisible par 10 selon Y
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'I0',[I0]) 
      
    d6run(d6OutFolder,newConfigFile)



#%%

#%
# Load the contents from the file, which creates a new dictionary with all experimental infos
in_file = open(JSONpath,"r")
dictExp = json.load(in_file)   


#generate the list of  files
lExp=[]
for d in dictExp:
    lExp.append(d)        
lExp=np.sort(lExp) 

import time
t=time.time()

w= 0.06
#All constant avalanche friction runs for figure 8 and 10

for j in [0, 1, 2, 3, 4, 5, 6, 7, 8]:   
    sE=lExp[j] #Selected exeperiment
    sdictE=dictExp[sE]
    delta_y=w/6
    if j<4:
        muR, mu = 0.3, 0.75
    else:
        muR, mu = 0.23, 0.44
    
    if j==8:
        nsamples, substeps = 8, 120 
    else:
        nsamples, substeps = 6, 80
    rund6py(sdictE, delta_mu=0., muRigid = muR, mu=mu, prop='', fps=15, nFrames=int(sdictE['nFrames']//10+5), nSamples=nsamples, I0_start=0.005, delta_mu_start=0., rand=1, substeps=substeps, W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y, mudoor= muR)





# High quality simulations

# fps=150
# for j in [5,6]:
#     sE=lExp[j] #Selected exeperiment
#     sdictE=dictExp[sE]    
#     for w in [0.06]:
#         delta_y=w/6
#         if w>=0.1:
#             delta_y=0.01
#         if j<4:
#             muR, mu = 0.3, 0.75
#         else:
#             muR, mu = 0.23, 0.44
    
#         if j==8:
#             nsamples, substeps = 8, 120 
#         else:
#             nsamples, substeps = 6, 80
            
#         rund6py(sdictE, delta_mu=0., muRigid = muR, mu=0.44, prop=fignames[j]+'infNorm_fps150_tst', fps=fps, nFrames=int(sdictE['nFrames']//(fps/150)+5*(fps/15)), nSamples=nsamples, I0_start=0.005, delta_mu_start=0., rand=1, substeps=int(substeps/(fps/15)), W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y, mudoor= muR)
        
        # rund6py(sdictE, delta_mu=0.26, muRigid = muR, mu=0.38, prop=fignames[j]+'muI', fps=15, nFrames=int(sdictE['nFrames']//10+5), nSamples=nsamples, I0_start=0.005, delta_mu_start=0., rand=1, substeps=substeps, W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y, mudoor= muR)
        
#         rund6py(sdictE, delta_mu=0.0, muRigid = muR, mu=mu, prop=fignames[j]+'hyst_th', fps=15, nFrames=int(sdictE['nFrames']//10+5), nSamples=nsamples, I0_start=0.005, delta_mu_start=0.10, rand=1, substeps=substeps, W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y, mudoor= muR)

#         rund6py(sdictE, delta_mu=0.0, muRigid = muR, mu=mu, prop=fignames[j]+'hyst_th', fps=15, nFrames=int(sdictE['nFrames']//10+5), nSamples=nsamples, I0_start=0.005, delta_mu_start=0.15, rand=1, substeps=substeps, W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y, mudoor= muR)

#         rund6py(sdictE, delta_mu=0.0, muRigid = muR, mu=mu, prop=fignames[j]+'hyst_th', fps=15, nFrames=int(sdictE['nFrames']//10+5), nSamples=nsamples, I0_start=0.005, delta_mu_start=0.2, rand=1, substeps=substeps, W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y, mudoor= muR) 
           
# j=7
# sE=lExp[j] #Selected exeperiment
# sdictE=dictExp[sE]
# for w in [ 0.01, 0.02, 0.04, 0.1, 0.15, 0.2]:
#     delta_y=w/4
#     if w>=0.1:
#         delta_y=0.02
#     if w<=0.02:
#         delta_y=w/4  
#     muR, mu = 0.23, 0.44  
#     nsamples, substeps = 6, 80
#     rund6py(sdictE, delta_mu=0.0, muRigid = muR, mu=mu, prop=fignames[j]+'width_var', fps=15, nFrames=int(sdictE['nFrames']//10+5), nSamples=nsamples, I0_start=0.005, delta_mu_start=0.00, rand=1, substeps=substeps, W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y, mudoor= muR) 
    
# high friction on door
# j=7  
# sE=lExp[j] #Selected exeperiment
# sdictE=dictExp[sE]
# for w in [0.06]:
#     delta_y=w/6
#     if j<4:
#         muR, mu = 0.3, 0.75
#     else:
#         muR, mu = 0.18, 0.44
#     rund6py(sdictE, delta_mu=0., muRigid = muR, mu=mu, prop=fignames[j]+'mu_door_high', fps=15, nFrames=int(sdictE['nFrames']//10+5), nSamples=6, I0_start=0.005, delta_mu_start=0., rand=1, substeps=80, W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y, mudoor= 1)  


    
    

