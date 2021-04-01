
# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Author : Gauthier Rousseau
import opyf  # from opyflow library some rendering function may be employed
sys.path.append(
    '/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python/imageProcessing')
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
# print(r'\includegraphics{test_2.pdf}')
driveFolder = '/media/gauthier/Data-Gauthier/Gauthier'
maind6OutFolder = '/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/'
paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(
    maind6OutFolder, 4)
mainExpFolder = driveFolder + \
    '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs_opyf'
folderOut='/media/gauthier/Data-Gauthier/Gauthier/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular-all/dry-granular.wiki/collapses/fit_runout/'
plt.rcParams['text.usetex']= True

# %% load all the exp profile and detect the run out
#%% gravels
plt.close('all')
fig, ax = plt.subplots(1, 1)
xfC = []
H=[]
for Nrun in range(0,4):
    runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
    D=runExp1.expD[-1]
    X = runExp1.vecXexpD
    xf = X[np.where(D < 0.001)[0][0]]
    ax.plot(X, D)
    ax.plot(xf,0,'+')
    xfC.append(xf)
    H.append(runExp1.dictE['H'])
    
xfC = np.array(xfC)
H=np.array(H)

plt.close('all')
fig, ax = plt.subplots(1, 1)
for mu in np.linspace(0.45,0.75,7):
    angles = np.array([0, 5, 10,15])*np.pi/180
    b= 1/(mu-np.tan(angles))


    ax.plot(b, xfC / H, '-+',linewidth=1.5, label='$\mu=$' + str(np.round(mu,2)))
    p, cov = np.polyfit(b, xfC/H, 1, cov=True)
    err = np.sqrt(np.diag(cov))
    ax.plot(b,np.polyval(p,b),'--r',linewidth=1, label='fit $\mu=$' + str(np.round(mu,2)) + '- err='+toS(err[0]/p[0]*100,2)+' $\%$')
ax.legend()
ax.set_ylabel(r'$x_f/h_0$',fontsize=12)
ax.set_xlabel(r'$\frac{1}{\mu-tan(\theta)}$',fontsize=12)
# ax.plot(b,b)
ax.set_title(r'gravels - $\tan \theta_{start}$= 0.65')
fig.savefig(folderOut+'gravels-fit.svg')
#%% beads
plt.close('all')
fig, ax = plt.subplots(1, 1)
xfC = []
H=[]
for Nrun in range(4,7):
    runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
    D=runExp1.expD[-1]
    X = runExp1.vecXexpD
    xf = X[np.where(D < 0.001)[0][0]]
    
    ax.plot(X, D)
    ax.plot(xf,0,'+')
    xfC.append(xf)
    H.append(runExp1.dictE['H'])
xfC = np.array(xfC)
H=np.array(H)

fig, ax = plt.subplots(1, 1)
xfC = np.append(xfC, 0.55)
H=np.append(H,H[-1])
for mu in np.linspace(0.35,0.55,5):
    angles = np.array([0, 5, 10,15])*np.pi/180
    b= 1/(mu-np.tan(angles))

    

    ax.plot(b, xfC / H, '-+',linewidth=1.5, label='$\mu=$' + str(np.round(mu,2)))
    p, cov = np.polyfit(b, xfC/H, 1, cov=True)
    err = np.sqrt(np.diag(cov))
    ax.plot(b,np.polyval(p,b),'--r',linewidth=1, label='fit $\mu=$' + str(np.round(mu,2)) + '- err='+toS(err[0]/p[0]*100,2)+' $\%$')

ax.legend()
ax.set_ylabel(r'$x_f/h_0$',fontsize=12)
ax.set_xlabel(r'$\frac{1}{\mu-tan(\theta)}$', fontsize=12)
ax.set_title(r'beads- $\tan \theta_{start}$= 0.43 ')

fig.savefig(folderOut + 'beads-fit.svg')

#%% compare with simu 

R1 = d6py.NumericalRun('/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/Run_00_3D_Door_mu=0.7_muRigid=0.0_W_10.0cm_gravels-2.7mm_Slope=0deg_delta_mu=0.00_substeps_40_fracH=0.8_I0_start=0.0000_delta_mu_start=0.00test-runout')
R2 = d6py.NumericalRun('/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/Run_01_3D_Door_mu=0.7_muRigid=0.0_W_10.0cm_gravels-2.7mm_Slope=5deg_delta_mu=0.00_substeps_40_fracH=0.8_I0_start=0.0000_delta_mu_start=0.00test-runout/')
R3 = d6py.NumericalRun('/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/Run_02_3D_Door_mu=0.7_muRigid=0.0_W_10.0cm_gravels-2.7mm_Slope=10deg_delta_mu=0.00_substeps_40_fracH=0.8_I0_start=0.0000_delta_mu_start=0.00test-runout/')
R4 = d6py.NumericalRun('/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/Run_03_3D_Door_mu=0.7_muRigid=0.0_W_10.0cm_gravels-2.7mm_Slope=15deg_delta_mu=0.00_substeps_40_fracH=0.8_I0_start=0.0000_delta_mu_start=0.00test-runout/')
selectedRuns=[R1, R2, R3, R4]
# selectedRuns = [Ref[0], selectedRuns[0]]
plt.close('all')
fig, ax = plt.subplots(1, 1)
xfM = []
HM=[]
for sR in selectedRuns:
    nF = int(sR.dConfig['nFrames'])

    sR.loadVTK(nF)


# sR.plotContour(ax)
    sR.h=np.zeros(len (sR.grid_x[:,sR.nYplot,0]))
    for i in range(len(sR.grid_x[:, sR.nYplot, 0])):
        ind = np.where(sR.reshaped_Phi[i, sR.nYplot,:] > 0.5)[0]
        if len(ind) > 0:        
            sR.h[i] = sR.grid_z[0, 0, ind[-1]]
        else:
            sR.h[i] = np.nan
    X = sR.grid_x[:, sR.nYplot, 0]
    xfM = X[np.where(sR.h < 0.001)[0][0]]
    HM.append(0.8 * sR.dConfig['box'][2])
xf = np.array(xfM)
H = np.array(H)

for mu in np.linspace(0.45,0.75,10):
    angles = np.array([0, 5, 10,15])*np.pi/180
    b= 1/(mu-np.tan(angles))

    

    ax.plot(b, xfC / H, '-+',linewidth=1.5, label='$\mu=$' + str(np.round(mu,2)))
    p, cov = np.polyfit(b, xfC/H, 1, cov=True)
    err = np.sqrt(np.diag(cov))
    ax.plot(b, np.polyval(p, b), '--r', linewidth=1, label='fit $\mu=$' + str(np.round(mu, 2)) + '- err=' + toS(err[0] / p[0] * 100, 2) + ' $\%$')

ax.legend()
ax.set_ylabel(r'$x_f/h_0$',fontsize=12)
ax.set_xlabel(r'$\frac{1}{\mu-tan(\theta)}$', fontsize=12)
ax.set_title(r'gravels- $\tan \theta_{start}$= 0.65 - from model to fit')
fig.savefig(folderOut + 'gravels-fit-from-model.svg')
#%%



R1 = d6py.NumericalRun('/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/Run_04_3D_Door_mu=0.38_muRigid=0.0_W_8.0cm_glass-beads-0.47mm_Slope=0deg_delta_mu=0.00_substeps_40_fracH=0.8_I0_start=0.0000_delta_mu_start=0.00talush/')

#%%
plt.close('all')
fig, ax = plt.subplots(1, 1)
R1.loadVTK(8)
R1.opyfPointCloudColoredScatter(ax)
nF = int(R1.dConfig['nFrames'])
masks_low = []
for i in range(5, 6):
    R1.loadVTK(i)
    R1.opyfPointCloudColoredScatter(ax, vmin=0, vmax=0.01)
ax.set_aspect('equal', adjustable='box') 
plt.show()
#%%
    R1.calculateNormVelocity()
    R1.IMvel = np.flipud(R1.normV.T)
    mask=np.zeros_like(R1.IMvel[:,4,:])
    ax.cla()
    ax.set_aspect('equal', adjustable='box')

    mask[np.where(R1.IMvel[:,4,:]> 0.008)] = 1

    masks_low.append(mask)
mask_final  = np.zeros_like(mask)
# fig.axes[-1].remove()
for m in masks_low:
    mask_final[np.where(m == 1.)] = 1
ax.set_aspect('equal', adjustable='box')    
ax.imshow(mask_final,cmap='hot',interpolation='gaussian',extent=R1.extentR)
fig.show()

#%%