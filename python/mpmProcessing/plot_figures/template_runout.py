#%%
import sys
import opyf  # from opyflow library some rendering function may be employed
sys.path.append(
    './..')
sys.path.append(
    './../..')
import d6py
import sys
sys.path.append(
    './../../imageProcessing')

from Tools_collapses import mask_collapses, mask_collapses2
import matplotlib.pyplot as plt
from d6py.Tools import *
import sys
import numpy as np

plt.rcParams['font.size'] = 8.0
plt.rcParams['xtick.labelsize'] = 7.0
plt.rcParams['ytick.labelsize'] = 7.0
plt.rcParams['ytick.labelsize'] = 8.0
plt.rcParams['axes.linewidth'] = 0.8


plt.close('all')

def generate_fig():

    fig, ax = plt.subplots(1,1)

    fig.set_size_inches(3.,2.2)
    ax2=fig.add_axes([0.1,0.9,0.8,0.1])
    ax.set_position([0.13, 0.15, 0.83, 0.72])
    ax.set_xlabel('$t$ [s]')
    ax2.grid()
    ax2.set_axis_off()
    ax2.plot([10,10],[10,11], ':',lw=0.9,c='k',label='Experience')
    ax2.plot([10,10],[10,11], '-',lw=0.9,c='k',label='Simulations - $\mu =0.44$')
    ax2.set_ylim([0,1])
    ax2.set_xlim([0,1])
    ax2.legend(edgecolor='w',facecolor='w',ncol=2,bbox_to_anchor=(0.7, 0.5, 0.4, 0.2),fontsize=leg_fontsize, framealpha=0)
    
    return fig, ax, ax2
#%%
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth   
#%%


def extractRunOut_time(sR,iFrame,iframe_deb=3,th=0.001):
    sR.loadVTK(iFrame)
    fps=sR.dConfig['fps']
    sR.h=np.zeros(len (sR.grid_x[:,sR.nYplot,0]))
    tf=(iFrame-iframe_deb)/fps
    contours=d6py.Tools.findContours(sR.grid_x[:,0, 0], sR.grid_z[0,0,:], sR.reshaped_Phi[:,sR.nYplot,:], 0.5)
    X=sR.grid_x[:, sR.nYplot, 0]
    for i in range(len(X)):
        ind=np.where(X[i]>contours[0][:,0])[0]
        if len(ind) > 0:        
            sR.h[i] = contours[0][ind[-1],1]
        else:
            sR.h[i] = np.nan
    ind0=np.where(sR.h > th)[0]
    if len(ind0)>0:  
        xfM = X[ind0[-1]]
    else:
        xfM = X[-1]
    
    return xfM, tf
    
    
    
def loadRunOutMod(sR,vth=0.01,hth=0.005):
    xfMss, VxfMss, tfs, maxVss = [], [], [], []
    nF = int(sR.dConfig['nFrames'])
    fps=sR.dConfig['fps']
    sR.loadVTK(3)
    k = 0
    sR.nYplot=2
    sR.h=np.zeros(len (sR.grid_x[:,sR.nYplot,0]))
    xfMs, maxVs, times = [], [], []
    for i in range(3,nF):
        sR.loadVTK(i)
        sR.calculateNormVelocity()
        indS=np.where(sR.reshaped_Phi[:,sR.nYplot,:]>0.99)
        norm=sR.normV[:,sR.nYplot,:]
        ind_HV=np.where(norm[indS]>vth)
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
        
        # xfM = X[np.where(sR.h < hth)[0][0]]
        sel_pts_higher_than_bed=sR.pointsp[np.where(sR.pointsp[:,2]>hth),0]
        sorted_index_array = np.argsort(sel_pts_higher_than_bed[:])
        xfM=np.mean(sel_pts_higher_than_bed[0,sorted_index_array[0][-100:]])
        if len(xfMs)>0:
            if xfM>=xfMs[-1]:
                xfMs.append(xfM)
            else: 
                xfMs.append(xfMs[-1])
        else:
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
    
    return xfMss, VxfMss, tfs, maxVss, times, dt
#%%  
    
def  loadRunOutExp(runExp,hth=0.005,vth=0.01):  
    mainExpFolder = driveFolder + \
            '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/' 
    listExpFolders=np.sort([mainExpFolder+'/'+d for d in runExp.dictExp ])
    Nrun=int(runExp.runNumber)
    Time, PointsSelected, VelSelected=opyf.hdf5_ReadUnstructured2DTimeserie(listExpFolders[Nrun]+'/Run_'+format(Nrun,'02.0f')+'_unstructured_good_features_to_track_velocities.hdf5')
    shiftExp=5

    #%
    k = 0
    plt.close('all')
    X = runExp.vecXexpD
    xfEs, timesE, maxVsE= [], [], []
    nf=len(runExp.expD)+ 30 + shiftExp-1
    for i in np.arange(30+shiftExp,nf):
        D=runExp.expD[i - 30 - shiftExp]
        ind = np.where(D < hth)
        if len(ind[0])<1:
            xf=np.nan
        else:
            xf = X[ind[0][0]]
        norm=(runExp.Ux[i - 30 - shiftExp]**2+runExp.Uy[i - 30 - shiftExp]**2)**0.5
        # maxVsE.append(np.nanmax(norm))
        maxVsE.append(np.nanmax(VelSelected[i - 30 - shiftExp]))
        if len(np.where(norm>vth)[0])>0:
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
    xfEs = smooth(xfEs,3)
    maxVsE = np.array(maxVsE)
    maxVsE = smooth(maxVsE,3)

    timesE = np.array(timesE)
    dtE=timesE[1]-timesE[0]
    VxfEs=(xfEs[1:]-xfEs[:-1])/dtE 
    
    return xfEs, timesE, maxVsE, VxfEs, tf_exp, timesE, dtE
# %%


