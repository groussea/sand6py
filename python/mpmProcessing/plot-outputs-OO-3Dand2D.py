#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:52:28 2019

@author: Gauthier Rousseau
"""


import sys
import numpy as np
import matplotlib.pyplot as plt
driveFolder='/scratch/garousse/'
driveFolder='/media/gauthier/Gauthier_Backup/'
#driveFolder='/media/garousse/Gauthier_Backup/'
sys.path.append(driveFolder+'TAF/TAF_EPFL/current_work/OPyF-Project/github/opyFlow/')
d6Path='/media/gauthier/Data-Gauthier/programs/gitLab/sand6/build'
# d6Path='/scratch/garousse/TAF/TAF_inria/INRIA_current_work/GitLab/sand6/build'
sys.path.append(d6Path+'/../python')
import opyf
import d6py

maind6OutFolder=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs/Tests'

paths,folders,listDictConf,listNumRun=d6py.findOutSand6Paths(maind6OutFolder,4)


#%%
N=7

selectedRuns0,selectedDict=d6py.whereSand6OutFromParms(listNumRun,runNumber=N,dimSim=3,delta_mu_start=0.0,delta_mu=0)
selectedRuns1,selectedDict=d6py.whereSand6OutFromParms(listNumRun,runNumber=N,dimSim=3,I0_start=0.004,delta_mu_start=0.07,delta_mu=0)

# selectedRuns2,selectedDict=d6py.whereSand6OutFromParms(listNumRun,runNumber=N,dimSim=3,delta_mu_start=0.09)
# selectedRuns2,selectedDict=d6py.whereSand6OutFromParms(listNumRun,runNumber=N,dimSim=3,delta_mu_start=0.09,delta_mu=0,I0_start=0.016)


selectedRuns=[selectedRuns0[-4],selectedRuns0[-3],selectedRuns0[-2]]
# selectedRuns=selectedRuns1[-1:]

for sR in selectedRuns:
    sR.scLength(0.01)
mainExpFolder=driveFolder+'/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs_opyf'
runExp1=d6py.ExperimentalRun(N,mainExpFolder,loadField=True)

runExp1.scLength(0.01)


plt.rc('text', usetex=True)#from opyflow library some rendering function may be employed

vec_cap_time=np.array([0,1,2,4,8,16,32,64,128])*runExp1.typicalTime


mainSaveFolder=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs/Tests/outputs/'
listSaveFolders=np.sort([mainSaveFolder+d for d in runExp1.dictExp])
for l in listSaveFolders:
    opyf.mkdir2(l)
cmap=opyf.custom_cmap.make_cmap_customized(Palette='green')
#%%
plt.close('all')
fig, [ax,ax1] = plt.subplots(2,1)
#d6Out_title= ''.join([(s!='_')*s for s in d6Out])
#ax.set_title(d6Out_title,fontsize=9)
kt=0
plt.ion()
for ifile in range(10,220,10):
    fig.set_size_inches((8*2,4.5*2))
    ax.cla()
    ax1.cla()
    h1, l1 = [], []
    ls=['--','-','-.','-']
    c=[0.2,0.8,1,0.1]
    for sR,i in zip(selectedRuns,range(len(selectedRuns))):
        sR.loadVTK(int(ifile*sR.dConfig['fps']//150))
        # sR.plotVelocity(ax,vmin=0,vmax=.5)
        # sR.plotI(ax,vmin=0,vmax=.1)
#        sR.plotPressure(ax,vmin=0,vmax=1000)
        # ax.imshow(np.mean(sR.P[:,:,:],axis=1),extent=sR.extentR,vmin=0,vmax=2000)
        # ax.imshow(np.mean(sR.epsilon21[:,:,:],axis=1),extent=sR.extentR,vmin=-50,vmax=50)
        # sR.plotContourVel(ax,levels=[0.1],linewidths=1,linestyles=ls[i%4],colors=[cmap((i+1)/(len(selectedRuns)+1))]) 
        sR.plotContour(ax,levels=[0.5],linewidths=1,linestyles=ls[i%4],colors=[cmap((i+1)/(len(selectedRuns)+1))]) 
  
        if i==0:
            # sR.plotI(ax,vmin=0,vmax=.1)
            sR.opyfPointCloudColoredScatter(ax,vmin=0,vmax=.5,s=4,cmap=cmap)
        h,l = sR.CS.legend_elements(str(sR.dimSim)+"D- \mu_{rb}="+format(sR.dConfig['muRigid'],'0.2f') +"- frac_h="+format(sR.dConfig['frac_h'],'0.2f') + "- \mu= "+ format(sR.dConfig['mu'],'0.2f')+ "- I0-start= "+ format(sR.dConfig['I0_start']) +"- \delta \mu_s= "+ format(sR.dConfig['delta_mu_start'],'0.4f')+" - \phi")        
        h1, l1=h1+h, l1+l
    runExp1.plotField(ax1,(ifile-7),vmin=0,vmax=.5,cmap=cmap)
    runExp1.plotDepthProfile(ax1,ifile-7,linestyle='--',color='k',linewidth=1,label="Experience")
    runExp1.plotDepthProfile(ax,ifile-7,linestyle='--',color='k',linewidth=1,label="Experience")
    ax.legend(h1+runExp1.l, l1+[r"Expe - tan($\theta$)=0.43 $\pm$ 0.03"],fontsize=8)
    ax.axis('equal')
    ax1.axis('equal')
    ax.set_ylim([-2,12])
    ax.set_xlim([-22.6,50])
    ax1.set_ylim([-2,12])
    ax1.set_xlim([-22.6,50])


#    
#    d6py.setFigure(fig,ax1,runExp1.scaleLength,unit='cm')
#    d6py.setFigure(fig,ax1,runExp1.scaleLength,unit='cm')
#    ax.set_ylim(np.array([-0.01,runExp1.dictE['H']+0.01])/runExp1.scaleLength)
#    ax.set_xlim(np.array([-runExp1.dictE['L'],runExp1.xmax*runExp1.dictE['H']])/runExp1.scaleLength)
#    ax.text(20,6.5,s='Time='+format(runExp1.Time[ifile],'10.2f')+' s', fontsize=8)

    plt.show()
    plt.pause(0.5)
#
#    if runExp1.Time[ifile]>=vec_cap_time[kt]:
#        kt+=1
    outFolder=listSaveFolders[runExp1.runNumber]+'/'+sR.dConfig['folder']
    opyf.mkdir2(outFolder)
    fig.savefig(outFolder+'/test2_compMPM_EXP_t='+format(runExp1.Time[ifile],'02.2f')+'.png')
#%%
#import subprocess
#import os
#for p in paths:
#    os.chdir(p)
#    cmd = str('rm -r frame*')
#    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
#    (out, err) = proc.communicate()
#    print("program output:", out)
#
#    

