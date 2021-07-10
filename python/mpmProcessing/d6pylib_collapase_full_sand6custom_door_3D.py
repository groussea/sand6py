
#!/usr/bin/env python
"""
Created on Mon Oct 16 14:35:06 2017


Python script to run the several cases of JFM
readconfigFile(configFile) return a dictionary 

@author: Gauthier
"""

import numpy as np  
import json
import sys, os


fileDir = os.path.dirname(os.path.abspath(__file__))

print(fileDir)

d6Path=fileDir+'/../../build'



mainOutFolder=d6Path+'/out'



JSONpath=fileDir+'/../Granular_Collapses_Experimental_Informations.json'

os.chdir(d6Path)
#add d6py module pythonpath
sys.path.append(d6Path+'/../python')



import d6py
from d6py.d6python3D import * # python must be reload if d6python2D was imported

#% by default, mainOutPut folder is in the build folder but it is highly recommended to set your own mainOutPut folder since outputs are generally large
# mainOutFolder=d6Path+'/out'
# mainOutFolder='/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/'
mainOutFolder='/home/gauthier/sorties_sand6/'
d6py.mkdir2(mainOutFolder) 

#%%
def rund6py(sdictE,**args):

    fracH=0.8
    
    if sdictE['camType']=='Phantom':
        sdictE['L']=sdictE['L']+0.1
        sdictE['Ltot']=sdictE['Ltot']+0.1
    L=sdictE['L']+0.01    
    Lmod=sdictE['Ltot']
    nFrames=args.get('nFrames',int(sdictE['nFrames']/10))
    if (sdictE['camType']=='BW') & (sdictE['Slope']==15.):
        Lmod = sdictE['Ltot'] + 0.5
        
    if (sdictE['camType']=='BW') & (sdictE['Slope']==20.):
        Lmod = sdictE['Ltot'] + 0.8
        
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
    prop=args.get('prop','test')
    rand = args.get('rand', 1)
    W = args.get('W', 0.08)

    if door=='with':
        runName=str('Run_'+format(j,'02.0f')+'_3D_Door_mu='+str(mu)+'_muRigid='+str(muRigid)+'_W_'+format(W*100,'.1f')+'cm_'+sdictE['grainType']+'_Slope='+format(sdictE['Slope'],'.0f')+'deg_delta_mu='+format(delta_mu,'.2f')+'_substeps_'+str(substeps)+'_fracH='+str(fracH)+'_I0_start='+format(I0_start,'.4f')+'_delta_mu_start='+format(delta_mu_start,'.2f')+prop)    
    else:
        runName=str('Run_'+format(j,'02.0f')+'_3D_no_Door_start_at_0.13_s_'+sdictE['grainType']+'_Slope='+format(sdictE['Slope'],'.0f')+'deg_delta_mu='+format(delta_mu,'.3f')+'_substeps_'+str(substeps)+'_fracH='+str(fracH)+'_I0_start='+format(I0_start,'.4f')+'_delta_mu_start='+format(delta_mu_start,'.2f')+prop)      
       

    d6OutFolder=mainOutFolder+'/'+runName
    d6py.mkdir2(d6OutFolder) 
    
    newConfigFile=d6OutFolder+'/collapse.3d_'+runName+'m.conf'
    configFilein=d6Path+'/../scenes/collapse.3d.LHE.Door.conf'

    d6py.modifyConfigFile(configFilein,newConfigFile,'box',[Lmod, W,Hmod])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'fps',[args.get('fps',15)])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'gravity',[+9.81*np.sin(sdictE['Slope']*np.pi/180), 0,-9.81*np.cos(sdictE['Slope']*np.pi/180)])
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
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.06 veldoor:0.8 ts:'+format(ts,'1.2f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L)+' wsw:'+str(wsw))
        else:
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.13 veldoor:0.7 ts:'+format(ts,'1.2f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L)+' wsw:'+str(wsw))
    else:
        d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.0001 veldoor:100 ts:'+format(ts,'1.2f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L))

    d6py.modifyConfigFile(newConfigFile,newConfigFile,'boundary','top:slip left:slip right:free front:slip back:slip bottom:stick')
    delta_x=args.get('delta_x',0.005)
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


# for j in range(8, 9):
for j in [8]:   
    sE=lExp[j] #Selected exeperiment
    sdictE=dictExp[sE]
    for w in [0.01,0.06,0.12]:
        # rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=np.round(sdictE['mu'] + 0.01, 2), prop='test', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=0, rand=0, substeps=40,W=w,wsw=0.01)
        delta_y=w/6
        rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=0.44, prop='resZ40_HRx', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=3, I0_start=0.005, delta_mu_start=0., rand=0, substeps=80, W=w+2*delta_y, wsw=delta_y, I0=0.0279,delta_y=delta_y)

for j in [7]:   
    sE=lExp[j] #Selected exeperiment
    sdictE=dictExp[sE]
    for w in [0.06]:     
        delta_y=w/6   
        rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=0.54, prop='hystResZ40', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=3, I0_start=0.005, delta_mu_start=0.1, rand=0, substeps=60, W=w+2*delta_y, wsw=delta_y, I0=0.0279,delta_y=delta_y)  
        rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=0.38, prop='38resZ40', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=3, I0_start=0.005, delta_mu_start=0., rand=0, substeps=60, W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y)    
        rund6py(sdictE, delta_mu=0.26, muRigid=0.18, mu=0.38, prop='muIresZ40', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=3, I0_start=0.005, delta_mu_start=0., rand=0, substeps=60, W=w+2*delta_y, wsw=delta_y, I0=0.279,delta_y=delta_y)       
    # for dmu,dmus in zip([0.05],[0.]):
    #     rund6py(sdictE, delta_mu=0.26, muRigid=0.18, mu=np.round(sdictE['mu'] + dmu, 2), prop='test-scaling', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=dmus, rand=1, substeps=80)

    # for dmu,dmus in zip([-0.05],[0.]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=np.round(sdictE['mu'] + dmu, 2), prop='test-scaling', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=dmus, rand=1, substeps=40)

    # for dmu,dmus in zip([0.15],[0.1]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.005, delta_mu_start=dmus, rand=1, substeps=40)

    # for dmu,dmus in zip([0.2],[0.15]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.005, delta_mu_start=dmus, rand=1, substeps=40)
    # for dmu,dmus in zip([0.2],[0.15]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.005, delta_mu_start=dmus, rand=1, substeps=40)


    # for dmu,dmus in zip([0.05],[0.]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-scaling', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=dmus, rand=1, substeps=40)

    # for dmu,dmus in zip([0.0],[0.]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-scaling', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=dmus, rand=1, substeps=40)

    # for dmu,dmus in zip([0.15],[0.1]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.005, delta_mu_start=dmus, rand=1, substeps=40)

    # for dmu,dmus in zip([0.1],[0.05]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.005, delta_mu_start=dmus, rand=1, substeps=40)

    # for dmu,dmus in zip([-0.05],[0.]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-scaling', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=dmus, rand=1, substeps=40)




### runs from 0 to 3 
    # for dmu,dmus in zip([0.1],[0.]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=dmus, rand=1, substeps=40)

    # for dmu,dmus in zip([0.],[0.]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=dmus, rand=1, substeps=40)

    # for dmu,dmus in zip([0.15],[0.15]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.005, delta_mu_start=dmus, rand=1, substeps=40)

    # for dmu,dmus in zip([0.0],[0.1]):
    #     rund6py(sdictE, delta_mu=0., muRigid=0., mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.005, delta_mu_start=dmus, rand=1, substeps=40)

#     for dmu,dmus in zip([0.],[0.]):
#         rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=dmus, rand=1, substeps=40)

#     for dmu,dmus in zip([-0.1],[0.]):
#         rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.00, delta_mu_start=dmus, rand=1, substeps=40)

#     for dmu,dmus in zip([0.05],[0.15]):
#         rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.005, delta_mu_start=dmus, rand=1, substeps=40)

#     for dmu,dmus in zip([0.0],[0.1]):
#         rund6py(sdictE, delta_mu=0., muRigid=0.18, mu=np.round(sdictE['mu'] + dmu, 2), prop='test-door', fps=15, nFrames=int(sdictE['nFrames'] / 10), nSamples=2, I0_start=0.005, delta_mu_start=dmus, rand=1, substeps=40)
# ###


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

    
    
    
    

