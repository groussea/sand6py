# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Author : Gauthier Rousseau
import matplotlib
import sys, os
fileDir = os.path.dirname(os.path.abspath(os.__file__))

get_ipython().run_line_magic('matplotlib', 'qt5')

import opyf  # from opyflow library some rendering function may be employed

sys.path.append('/media/gauthier/DataSSD/programs/gitLab/sand6/python/mpmProcessing/plot_figures')
from template_runout import *



sys.path.append(
    '/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python/imageProcessing')
sys.path.append(
    '/media/gauthier/DataSSD/programs/gitLab/sand6/python/imageProcessing')
sys.path.append(
    '/media/gauthier/DataSSD/programs/gitLab/sand6/python/')
        
import matplotlib.pyplot as plt
from d6py.Tools import *
import d6py

# intialize exteral packages
import sys
import os
import numpy as np
import matplotlib

plt.rcParams['font.size'] = 8.0
plt.rcParams['xtick.labelsize'] = 7.0
plt.rcParams['ytick.labelsize'] = 7.0
plt.rcParams['ytick.labelsize'] = 8.0
plt.rcParams['axes.linewidth'] = 0.8


print(fileDir)

d6Path=fileDir+'/../../build-2d'

mainOutFolder=d6Path+'/out'

JSONpath=fileDir+'/../Granular_Collapses_Experimental_Informations.json'



driveFolder = '/media/gauthier/Data-Gauthier/Gauthier'
driveFolder = '/media/gauthier/DataSSD'

maind6OutFolder = '/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/'
maind6OutFolder = '/media/gauthier/DataSSD/sand6_out/'

paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths( maind6OutFolder, 4)
mainExpFolder = driveFolder + \
    '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'
Nrun = 7
scale = 0.01  # 1cm
runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
runExp1.scLength(scale)
mu = runExp1.dictE['mu']

#%%
# Run7

## 3D runs

R_3d, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=3, delta_mu=0.0,delta_mu_start=0,keyWord='width_var')
R2, selectedDict = d6py.whereSand6OutFromParms(listNumRun, runNumber=Nrun, keyWord='W_8',mu=0.44, delta_mu_start=0)

# resorted
R_2d, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=2, delta_mu=0.0,delta_mu_start=0,keyWord='resZ30')

R_3d = [R_3d[0], R_3d[-2], R_3d[-1], R2[0], R_3d[2]]

R_2d = [R_2d[-2], R_2d[-1], R_2d[3], R_2d[0], R_2d[1], R_2d[2]]
#
#%%
# init Vini
ifile=3
Vini = np.zeros((len(R_3d)))
for sR, i in zip(R_3d, range(len(R_3d))):
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    sR.nYplot=int(sR.resY//2)
    contrs = sR.findContourPhi(level=0.5)
    V = area(contrs[0])
    Vini[i] = V
# %% init andload all final sates
for sR in R_3d:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])
    print(nF)
    sR.nYplot=int(sR.resY//2)
    sR.loadVTK(int(nF* sR.dConfig['fps'] / 15))


for sR in R_2d:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])
    sR.nYplot=int(sR.resY//2)
    sR.loadVTK(int(nF * sR.dConfig['fps'] / 15))


#%%
plt.close('all')
fig, ax = plt.subplots(1,1)
mus=[]
hs=[]
for sR in R_2d:
    hmax=[]
    for k in range(10):
        hmax.append(np.mean(sR.grid_y[0,np.where(sR.reshaped_Phi[k,:]>0.5)[0][-1]]))
    # sR.hmax=np.mean(hmax)
    sR.hmax=np.max(sR.pointsp[:,1])
    print(sR.hmax)
    ax.plot(sR.hmax,sR.dConfig['mu'],'k+')
    mus.append(sR.dConfig['mu'])
    hs.append(sR.hmax)
ax.set_ylabel(r'$\mu$')
ax.set_xlabel(r'$h_{max}$')
#%%
from scipy import optimize
def moinslogvraisemblance(x):
    logvr=[]
    for k in range(len(hs)):
        logvr.append(1/2*(np.log(np.sqrt(2*np.pi))+np.log(1)+(((mus[k]-(x[0]+x[1]*hs[k]**x[2]))/(1))**2)))
    lvrai=np.nansum(logvr)
    return lvrai


# xopt = optimize.fmin(func=moinslogvraisemblance, x0=[1,1])
xopt = optimize.minimize(moinslogvraisemblance, [1,1,1],method='Nelder-Mead')
print(xopt)
hs=np.array(hs)
ax.plot(hs,xopt.x[0]+xopt.x[1]*hs**xopt.x[2])

plt.show()




#%%
plt.rcParams["text.usetex"]=True
plt.close('all')
fig, ax = plt.subplots(1,1,figsize=(4,3))
Ws=[1,2,4,6,15]
mu_eqs,hmaxs = [], []
i=0
for sR,w in zip(R_3d,Ws):
    hmax=[]
    H0 = sR.dConfig['box'][2]*0.8
    for k in range(2,8):
        hmax.append(sR.grid_z[0,0,np.where(sR.reshaped_Phi[k,3,:]>0.5)[0][-1]])
    # sR.hmax=np.mean(hmax)
    xfM, t_f = extractRunOut_time(sR,int(sR.dConfig['nFrames']),iframe_deb=3)
    sR.hmax=np.nanmax(sR.h)
    
    # sR.hmax=np.max(sR.pointsp[:,2])
    print(sR.hmax)
    hmaxs.append(sR.hmax)
    mu_eq=xopt.x[0]+xopt.x[1]*sR.hmax**xopt.x[2]
    # print(sR.hmax)
    mu_eqs.append(mu_eq)
    # ax.plot(w,sR.hmax,'k+')
    #evaluate the lost
    V = area(sR.findContourPhi(level=0.5)[0])

    lost = (Vini[i]-V)/Vini[i]*100
    # print(lost)
    i+=1

Ws=np.array(Ws)
mu_eq_err=xopt.x[1]*(sR.hmax+0.004)**xopt.x[2]-(xopt.x[1]*(sR.hmax-0.004)**xopt.x[2])

# ax.plot(Ws*0.01/sR.dConfig['box'][2]*0.8,mu_eqs,'3', color='darkblue',label=r'3d simulation  $\mu=0.44$',ms=10)

ax.plot((Ws*0.01/H0)**(-1),mu_eqs,'p', color='darkblue',label=r'3d B15 AR=' + format(sR.dConfig['box'][2]*0.8/sR.dConfig['column_length'],'1.1f') + r' $\mu=0.44$',ms=4)

ax.errorbar((Ws*0.01/H0)**(-1),mu_eqs, yerr=0.02, color='k',fmt='none', markersize=8,markeredgewidth=1, capsize=2, )

ax.set_ylabel(r'$\mu_{2D,eq}$')
ax.set_xlabel(r'$H_0/W$ ')



ws_Ls=np.linspace(0.04,0.30,100)
muIonescu = 0.38 + 0.18 * 0.05 / ws_Ls
ax.plot((ws_Ls/H0)**(-1),muIonescu,'--',c='darkred',label=' Ionescu et al. (2015) Correction',lw=1.4)

ax.set_position([0.15,0.15,0.8,0.8])
#%% save data
filename = "/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/walls/AR0.5_15deg/AR0.5_15deg_3D.csv"
slope= np.absolute(sR.slope)
variables = [['Ws',Ws],['H0',np.ones(len(mu_eqs))*H0],['slope',np.ones(len(mu_eqs))*slope],['Hf',hmaxs],['mu_2D_eq',mu_eqs]]
opyf.write_csvScalar(filename, variables)



#%%

#add the AR2 points

filename = "/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/walls/AR2_0deg/AR2_0deg_3D.csv"

header, data= opyf.read_csv(filename)
Ws = data[:,0]
mu_eqs = data[:,4]
H0 = data[0,1]
ax.plot((Ws*0.01/H0)**(-1),mu_eqs,'p', color='cornflowerblue',label=r'3d B00 AR2  $\mu=0.44$',ms=4)

ax.errorbar((Ws*0.01/H0)**(-1),mu_eqs, yerr=0.02, color='cornflowerblue',fmt='none', markersize=8,markeredgewidth=1, capsize=2, )


ax.legend()
plt.show()

fig.savefig(driveFolder+"/programs/gitLab/dry-granular/doc/article/figures/mu2Deq_W.pdf", dpi=150)

#%%

def moinslogvraisemblance2(x):
    logvr=[]
    mus=np.array(
[0.5659759452221533,
 0.4710612618560802,
 0.45433547394426244,
 0.44209445338970355,
 0.43807308784497123])-0.44
    w=np.array([1,2,6,10,20])
    for k in range(len(w)):
        logvr.append((mus[k]-(x[0]/w[k]**x[1]))**2)
    logvr=np.array(logvr)
    lvrai=np.nansum(logvr)
    print(lvrai)
    return lvrai


# xopt2 = optimize.fmin(func=moinslogvraisemblance2, x0=[1,1])
xopt2 = optimize.minimize(moinslogvraisemblance2, [0.1,5],options={'gtol': 1e-6, 'disp': True})
print(xopt2)
Ws=np.array(Ws)
mus=np.array(mu_eqs)
# ax.errorbar(Ws,mus, yerr=np.ones(len(Hs)), color='k',fmt='none', markersize=8,markeredgewidth=1, capsize=10)



wl=np.linspace(1,20,100)
ax.plot(wl,0.44+xopt2.x[0]/wl**xopt2.x[1],'k--')
X,Y=wl[10],0.44+xopt2.x[0]/wl[10]**xopt2.x[1]

ax.annotate(r'$\mu_{2D,eq}(W)=0.44+0.13/W^{1.87}$',(X,Y),xytext=(X+4,Y+0.05),arrowprops=dict(color='k',arrowstyle="->", connectionstyle="arc3"),backgroundcolor='w',fontsize=11)
ax.set_ylabel(r'$\mu_{2D,eq}$')
ax.set_xlabel(r'$W~\mathrm{[cm]}$ ')

[x, y, X, Y] = ax.get_position().bounds
ax.set_position([0.15, 0.15, 0.8, 0.8])
plt.show()
fig.set_size_inches(4., 3)
folderOut='/media/gauthier/Data-Gauthier/Gauthier/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular-all/dry-granular/doc/article/images/used_images/'
fig.savefig(folderOut + 'mu_w_run7_15.pdf')

fig.savefig(folderOut + 'mu_w_run7_15.png')
#%%,'


#%% 
# 3d

for sR in R_3d:
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))


#%%



# for sR, i in zip(SR, range(len(SR))):

#     ax = ax1
#     sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
#     sR.nYplot=int(sR.resY//2)
#     contrs = sR.findContourPhi(level=0.5)
#     V = area(contrs[0])
#     Vini[i] = V
