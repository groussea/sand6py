
#%%
import opyf  # from opyflow library some rendering function may be employed
sys.path.append(
    '/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python/imageProcessing')
sys.path.append(
    '/media/gauthier/DataSSD/programs/gitLab/sand6/python/imageProcessing')
sys.path.append('/media/gauthier/DataSSD/programs/gitLab/sand6/python')
from Tools_collapses import mask_collapses, mask_collapses2
import matplotlib.pyplot as plt
from d6py.Tools import *
import d6py
# intialize exteral packages
import sys
import os
import numpy as np
import matplotlib
# sys.path.append('/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python')
# sys.path.append('/media/gauthier/Data-Gauthier/programs/gitHub/opyflow')
%matplotlib qt5
# sys.stdout = open(os.devnull, 'w')
# intialize font type and size
plt.rcParams['font.size'] = 8.0
plt.rcParams['xtick.labelsize'] = 7.0
plt.rcParams['ytick.labelsize'] = 7.0
plt.rcParams['ytick.labelsize'] = 8.0
plt.rcParams['axes.linewidth'] = 0.8
#%%

driveFolder = '/media/gauthier/Data-Gauthier/Gauthier'
driveFolder = '/media/gauthier/DataSSD'
maind6OutFolder = '/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/'
paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(
    maind6OutFolder, 4)

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth        

#%%
Nrun=2
run=str(Nrun)

fignames=['G00', 'G05', 'G10', 'G15', 'B00', 'B05', 'B10', 'B15', 'B20']
run=fignames[Nrun]
R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0., runNumber=Nrun, dimSim=3, delta_mu_start=0., mu=0.75, keyWord='try2')


R2, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0., runNumber=Nrun, dimSim=3, mu=0.75, viscosity=0.0)

# R3, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0., runNumber=Nrun, dimSim=3, delta_mu_start=0.1, mu=0.54, keyWord='visc')

# allRuns, sD = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0., dimSim=3, delta_mu_start=0., mu=0.75, keyWord='try2')

#%%
plt.close('all')
Runs=[R1[0]]
xfMss, VxfMss, tfs, maxVss = [], [], [], []
tf=0
for sR in Runs[:]:
    nF = int(sR.dConfig['nFrames'])
    fps=sR.dConfig['fps']
    sR.loadVTK(3)
    k = 0
    sR.h=np.zeros(len (sR.grid_x[:,sR.nYplot,0]))
    xfMs, maxVs, times = [], [], []
    for i in range(3,nF):
        sR.loadVTK(i)
        sR.calculateNormVelocity()
        indS=np.where(sR.reshaped_Phi[:,sR.nYplot,:]>0.99)
        norm=sR.normV[:,sR.nYplot,:]
        ind_HV=np.where(norm[indS]>0.01)
        maxV=np.max(norm[indS])
        # maxV=np.max(sR.vel)
        maxVs.append(maxV)
        if len(ind_HV[0])>1:
            tf=k/fps
        contours=d6py.Tools.findContours(sR.grid_x[:,0, 0], sR.grid_z[0,0,:], sR.reshaped_Phi[:,sR.nYplot,:], 0.5)
        X=sR.grid_x[:, sR.nYplot, 0]
        for i in range(len(X)):
            
            ind=np.where(X[i]>contours[0][:,0])[0]
            if len(ind) > 0:        
                sR.h[i] = contours[0][ind[-1],1]
            else:
                sR.h[i] = np.nan
        # for i in range(len(sR.grid_x[:, sR.nYplot, 0])):
        #     ind = np.where(sR.reshaped_Phi[i, sR.nYplot,:] > 0.5)[0]
        #     if len(ind) > 0:

        #         sR.h[i] = sR.grid_z[0, 0, ind[-1]]
        #     else:
        #         sR.h[i] = np.nan

        xfM = X[np.where(sR.h < 0.007)[0][0]]
        xfMs.append(xfM)
        times.append(k/fps)
        
        k+=1
    tfs.append(tf)
    xfMs=np.array(smooth(xfMs, 3))
    times=np.array(times)
    dt=times[1]-times[0]
    VxfMs=(xfMs[1:]-xfMs[:-1])/dt
    VxfMs[np.where(VxfMs<0)]=0
    xfMss.append(xfMs)
    VxfMss.append(VxfMs)
    maxVss.append(maxVs)
#%%
tfs_exp =[]

for Nrun in range(Nrun,Nrun+1):
    mainExpFolder = driveFolder + \
            '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'
    runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
    listExpFolders=np.sort([mainExpFolder+'/'+d for d in runExp1.dictExp ])
    Time, PointsSelected, VelSelected=opyf.hdf5_ReadUnstructured2DTimeserie(listExpFolders[Nrun]+'/Run_'+format(Nrun,'02.0f')+'_unstructured_good_features_to_track_velocities.hdf5')
    shiftExp=5

    #%
    k = 0
    plt.close('all')
    X = runExp1.vecXexpD
    xfEs, timesE, maxVsE = [], [], []
    nf=len(runExp1.expD)+ 30 + shiftExp-1
    for i in np.arange(30+shiftExp,nf):
        D=runExp1.expD[i - 30 - shiftExp]
        ind = np.where(D < 0.007)
        if len(ind[0])<1:
            xf=np.nan
        else:
            xf = X[ind[0][0]]
        norm=(runExp1.Ux[i - 30 - shiftExp]**2+runExp1.Uy[i - 30 - shiftExp]**2)**0.5
        # maxVsE.append(np.nanmax(norm))
        maxVsE.append(np.nanmax(VelSelected[i - 30 - shiftExp]))
        if len(np.where(norm>0.03)[0])>0:
            tf_exp=k/150
        # if xf<0:
        #     xf=np.nan
        # plt.clf()

        # plt.imshow(norm)
        # plt.pause(0.1)

        xfEs.append(xf)    
        timesE.append(k/150)
        k+=1
    
    xfEs = np.array(xfEs)
    xfEs = smooth(xfEs,5)
    maxVsE = np.array(maxVsE)
    maxVsE = smooth(maxVsE,5)

    timesE = np.array(timesE)
    dtE=timesE[1]-timesE[0]
    VxfEs=(xfEs[1:]-xfEs[:-1])/dtE 
    tfs_exp.append(tf_exp)


    
# VxfEs[np.where(VxfEs<0) ]=0
#%%


plt.close('all')
plt.figure()
plt.plot(times[:-2], xfMss[0][:-2],'+:',c='g',lw=0.9, label='Num. front position '+run+' norm')
# plt.plot(times[:-2], xfMss[1][:-2],'1--',c='y',lw=1.2, label='Num. front position '+run+' visc')
plt.plot(timesE[8:-10], xfEs[8:-10],'--',c='k',lw=1.4, label='Exp. front position '+run+' norm')
plt.plot(times[:-2:]+dt/2, VxfMss[0][:-1:],'+-',c='g',lw=0.9, label='Num. velocity '+run+' norm')
# plt.plot(times[:-2:]+dt/2, VxfMss[1][:-1:],'1-',c='y',lw=1.2, label='Num. velocity '+run+' visc')
plt.plot(timesE[10:-21:]+dtE/2, smooth(VxfEs[10:-20:],25),'-',c='k',lw=0.7, label='Exp. velocity '+run+'')
plt.show()

plt.legend(fontsize=8)
fig=plt.gcf()

fig.set_size_inches(5,4)

ax=plt.gca()
ax.set_position([0.1, 0.1, 0.85, 0.85])
ax.set_xlabel('time [s]')
ax.set_ylabel('position [m] -- velocitiy [m/s]')

fig.savefig('/media/gauthier/DataSSD/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular-all/dry-granular.wiki/20220512/'+run+'_front.png')

#%%

plt.close('all')
plt.figure()
plt.plot(times[:-2], maxVss[0][:-2],'+:',c='g',lw=0.9, label='Num. max vel. '+run+' norm')
# plt.plot(times[:-2], maxVss[1][:-2],'1--',c='y',lw=1.2, label='Num. max vel. '+run+' visc')
plt.plot(timesE[8:-10], maxVsE[8:-10],'--',c='k',lw=1.4, label='Exp. max vel. '+run+' norm')
plt.plot(times[:-2:]+dt/2, VxfMss[0][:-1:],'+-',c='g',lw=0.9, label='Num. velocity '+run+' norm')
# plt.plot(times[:-2:]+dt/2, VxfMss[1][:-1:],'1-',c='y',lw=1.2, label='Num. velocity '+run+' visc')
plt.plot(timesE[10:-21:]+dtE/2, smooth(VxfEs[10:-20:],25),'-',c='k',lw=0.7, label='Exp. velocity '+run+'')
plt.show()

plt.legend(fontsize=8)
fig=plt.gcf()

fig.set_size_inches(5,4)

ax=plt.gca()
ax.set_position([0.1, 0.1, 0.85, 0.85])
ax.set_xlabel('time [s]')
ax.set_ylabel('position [m] -- velocitiy [m/s]')

fig.savefig('/media/gauthier/DataSSD/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular-all/dry-granular.wiki/20220512/'+run+'_max_vel_visc.png')

# %%
