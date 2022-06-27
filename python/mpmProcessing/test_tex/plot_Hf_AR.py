# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Author : Gauthier Rousseau
from ctypes import alignment
import matplotlib
import sys, os
sys.path.append('/media/gauthier/DataSSD/programs/gitLab/sand6/python/mpmProcessing/plot_figures')
from template_runout import *
fileDir = os.path.dirname(os.path.abspath(os.__file__))

get_ipython().run_line_magic('matplotlib', 'qt5')

import opyf  # from opyflow library some rendering function may be employed


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
Nrun = 4
scale = 0.01  # 1cm
runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
runExp1.scLength(scale)
mu = runExp1.dictE['mu']

#%
# RunAR2

## 3D runs

# R_2d, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=2, delta_mu=0.0,delta_mu_start=0,keyWord='W_1.5')

# R_3d, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=3, delta_mu=0.0,delta_mu_start=0,keyWord='W_1.5_R', muRigid=0.00)
# del R_3d[1]


# R_3d, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=3, delta_mu=0.0,delta_mu_start=0,keyWord='W_24', muRigid=0.23)
# del R_3d[3]

# #
# R_3d=R_3d+R_3db

plt.rcParams["text.usetex"]=True

#%%
def drawPowerLaw(ax,p1,p2,shy,sl,inv=0,frmt='01.2f'):

    pm=np.mean([p1, p2])-p1
    dp=p2-p1
    points=np.array(([10**(p1),10**(p2)]) )       
    ax.plot(points,points**sl*10**(shy), '-r',lw=1)
    ax.plot(points,[points[1-inv]**sl*10**(shy), points[1-inv]**sl*10**(shy)], '-r',lw=1)
    ax.plot([points[inv], points[inv]],points**sl*10**(shy), '-r',lw=1)
    ax.text(points[0]*10**(pm),points[1-inv]**sl*10**(0.05*(1-2*inv))*10**(shy), '1' ,horizontalalignment='center',verticalalignment='center')
    ax.text(points[inv]*10**(-0.07*(1-2*inv)),points[0]**sl*10**(shy)*10**(sl*dp/2), format(sl,frmt),horizontalalignment='center',verticalalignment='center' )

plt.close('all')
fig, [ax, ax2] = plt.subplots(1,2,figsize=(5,2.5))
# fig2, ax2 = plt.subplots(1,1,figsize=(4,3))


R_3d, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=3, delta_mu=0.0,delta_mu_start=0,keyWord='no_fric_serie_30', muRigid=0.0)
ifile=3

Vini = np.zeros((len(R_3d)))
for sR, i in zip(R_3d, range(len(R_3d))):
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    sR.nYplot=int(sR.resY//2)
    contrs = sR.findContourPhi(level=0.1)
    V = area(contrs[0])
    Vini[i] = V

for sR in R_3d:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])
    # print(nF)
    sR.nYplot=int(sR.resY//2)
    sR.loadVTK(int(nF* sR.dConfig['fps'] / 15))
    
 
#%
# fig, [ax, ax2] = plt.subplots(1,2,figsize=(5,2.5))
mu_eqs=[]
hmaxs=[]
ARs=[]
H0s=[]
xfMs=[]
i=0
col='seagreen'
for sR in R_3d:
    H0=sR.dConfig['box'][2]*0.8
    hmax=[]
    L0=(sR.Orx-0.01)
    if H0/L0<20:
        xfM, t_f = extractRunOut_time(sR,int(sR.dConfig['nFrames']),iframe_deb=3,th=0.001)
        xfM=np.max(sR.pointsp[:,0])
        for k in range(2,4):
            hmax.append(sR.h[k])
        sR.hmax=np.nanmax(sR.h)
        # sR.hmax=np.max(sR.pointsp[:,2])
        ax.loglog(H0/L0,H0/sR.hmax,'2',c=col)
        #evaluate the lost
        V = area(sR.findContourPhi(level=0.1)[0])
    
        [marker2]=ax2.loglog(H0/L0,xfM/L0,'2',c=col)
        hmaxs.append(sR.hmax)
        ARs.append(H0/L0)
        H0s.append(H0)
        xfMs.append(xfM)
        # lost = (Vini[i]-V)/Vini[i]*100
        # print(lost)
        i+=1

hmaxs=np.array(hmaxs)
ARs=np.array(ARs)
H0s=np.array(H0s)
xfMs=np.array(xfMs)

inds1=np.where(ARs<7)
inds2=np.where(ARs<7)

fitW = np.polynomial.polynomial.polyfit(np.log(ARs[inds2]), np.log(H0s[inds2]/hmaxs[inds2]), [0, 1])

fitW2 = np.polynomial.polynomial.polyfit(np.log(ARs[inds1]), np.log(xfMs[inds1]/(H0s[inds1]/ARs[inds1])), [0, 1])

# print(fit)

# fit_fn = np.polynomial.Polynomial(fit)
ARs_fit = np.logspace(0, 1.3, 10)
[line3]=ax.plot(ARs_fit,np.exp(fitW[0])*ARs_fit**fitW[1], ':k',
        label=r"fit wide ($a < 7$): "+ format(np.exp(fitW[0]),'.1f') +r" $a^{"+ format(fitW[1],'.2f') +r"}$",lw=1)
# ax.plot(ARs_fit,np.exp(fit[0])*ARs_fit**0.6, '--r',lw=1)
[line4]=ax2.plot(ARs_fit[0:5],np.exp(fitW2[0])*ARs_fit[0:5]**fitW2[1], '--k',
        label=r"fit wide ($a < 7$): "+ format(np.exp(fitW2[0]),'.1f') +r" $a^{"+ format(fitW2[1],'.1f') +r"}$",lw=1)
inds1=np.where(ARs>6)
fitW2 = np.polynomial.polynomial.polyfit(np.log(ARs[inds1]), np.log(xfMs[inds1]/(H0s[inds1]/ARs[inds1])), [0, 1])

[line4]= ax2.plot(ARs_fit[5:10],np.exp(fitW2[0])*ARs_fit[5:10]**fitW2[1], ':k',
        label=r"fit wide ($a > 7$): "+ format(np.exp(fitW2[0]),'.1f') +r" $a^{"+ format(fitW2[1],'.2f') +r"}$",lw=2)

legends=[]

# ax2.plot(ARs_fit,np.exp(fit2[0])*ARs_fit**1, '--r',lw=1)
# ax2.plot(ARs_fit,np.exp(fit2[0])*ARs_fit**0.66, '--r',lw=1)

shy, sl, p1, p2= 0.4, 1, 0.2, 0.4

drawPowerLaw(ax2,p1,p2,shy,sl,frmt='1.0f',)

# 10**(sl*dp/2)
shy, sl, p1, p2=0.6, 0.65, 0.9, 1.1
drawPowerLaw(ax2,p1,p2,shy,sl,frmt='1.1f')

R_3d1, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=3, delta_mu=0.0,delta_mu_start=0,keyWord='HR_fric_serie_30', muRigid=0.23)



# del R_3d1[1]

Vini = np.zeros((len(R_3d1)))
for sR, i in zip(R_3d1, range(len(R_3d1))):
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    sR.nYplot=int(sR.resY//2)
    contrs = sR.findContourPhi(level=0.1)
    V = area(contrs[0])
    Vini[i] = V



# % init andload all final sates
for sR in R_3d1:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])
    # print(nF)
    sR.nYplot=int(sR.resY//2)
    sR.loadVTK(int(nF* sR.dConfig['fps'] / 15))



mu_eqs=[]
hmaxs=[]
ARs=[]
H0s=[]
xfMs=[]
i=0

col='firebrick'
for sR in R_3d1[:]:
    H0=sR.dConfig['box'][2]*0.8
    hmax=[]
    xfM, t_f = extractRunOut_time(sR,int(sR.dConfig['nFrames']),iframe_deb=3,th=0.001)
    xfM=np.max(sR.pointsp[:,0])
    for k in range(2,4):
        hmax.append(sR.h[k])
    sR.hmax=np.nanmax(sR.h)
    sR.hmax=np.max(sR.pointsp[:,2])
    L0=(sR.Orx-0.01)
    ax.loglog(H0/L0,H0/sR.hmax,'1',c=col)
    
    #evaluate the lost
    V = area(sR.findContourPhi(level=0.1)[0])
    
    [marker]=ax2.loglog(H0/L0,(xfM)/L0,'1',c=col)
    xfMs.append(xfM)
    hmaxs.append(sR.hmax)

    ARs.append(H0/L0)
    H0s.append(H0)
    lost = (Vini[i]-V)/Vini[i]*100
    print(lost)
    i+=1

hmaxs=np.array(hmaxs)
ARs=np.array(ARs)
H0s=np.array(H0s)
xfMs=np.array(xfMs)



inds1=np.where(ARs>1)

fit = np.polynomial.polynomial.polyfit(np.log(ARs[inds1]), np.log(H0s[inds1]/hmaxs[inds1]), [0, 1])

inds1=np.where(ARs>3)

fit2 = np.polynomial.polynomial.polyfit(np.log(ARs[inds1]), np.log(xfMs[inds1]/(H0s[inds1]/ARs[inds1])), [0, 1])

print(fit)

fit_fn = np.polynomial.Polynomial(fit)
ARs_fit = np.logspace(0, 1.5, 10)

[line1]=ax.plot(ARs_fit,np.exp(fit[0])*ARs_fit**fit[1], '-k',
        label=r"fit narrow: "+ format(np.exp(fit[0]),'.1f') +r" $a^{"+ format(fit[1],'.1f') +r"}$",lw=1)

[line2]=ax2.plot(ARs_fit,np.exp(fit2[0])*ARs_fit**fit2[1], '-k',
        label=r"fit narrow ($a>3$): "+ format(np.exp(fit2[0]),'.1f') +r" $a^{"+ format(fit2[1],'.1f') +r"}$",lw=1)
  
shy, sl, p1, p2= 0.02, 0.5, 0.4, 0.6

drawPowerLaw(ax,p1,p2,shy,sl,inv=1,frmt='1.1f')

shy, sl, p1, p2= 0.07, fit2[1], 0.7, 0.9

drawPowerLaw(ax2,p1,p2,shy,sl,inv=1,frmt='1.1f')


shy, sl, p1, p2= 0.2, 0.67, 0.2, 0.4

drawPowerLaw(ax,p1,p2,shy,sl,frmt='1.2f')



ax.axis('equal')

plt.show()

# R_3d, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=3, delta_mu=0.0,delta_mu_start=0,keyWord='W_24', muRigid=0.23)



ax.set_xlabel('$a$')
ax.set_ylabel('$H_0/H_{\infty}$')
ax.set_position([0.6,0.13,0.38,0.85])
ax.set_ylim(0.9,20)
ax.set_xlim(0.9,20)

ax2.set_xlabel('$a$')
ax2.set_ylabel('$(L_{\infty}-L_0)/L_0$')
ax2.set_ylim(0.9,30)
ax2.set_xlim(0.9,20)

ax2.set_position([0.1,0.13,0.38,0.85])

fig.set_size_inches(7,3)


fig.text(0.05,0.94,'(a)',fontsize=9)
fig.text(0.55,0.94,'(b)',fontsize=9)

l_new_1=[ r"wide - $W  \rightarrow \infty$", r"narrow - $W = 1$ cm",]

h_new_1=[marker2, marker]



      
fig.legend(h_new_1, l_new_1, fontsize=7, framealpha=0.5, bbox_to_anchor=(0.2, 0.47, 0.1, 0.5))
fig.legend(h_new_1, l_new_1, fontsize=7, framealpha=0.5, bbox_to_anchor=(0.47, 0.17, 0.5, 0.1))

# fig.legend([line1,line3],[leg1,leg3],  fontsize=7, framealpha=0.5, loc=2)

# ax2.legend(h_new_1, l_new_1, fontsize=7, framealpha=0.5, loc=2)


ax.legend(fontsize=7, loc=2)
ax2.legend(fontsize=7,loc=4)

plt.show()


fig.savefig(driveFolder+"/programs/gitLab/dry-granular/doc/article/figures/scalings_Hf_Lf_a.pdf", dpi=150)


#%%

fig2, ax3 = plt.subplots(1,1,figsize=(5,2.5))

ax3.loglog(ARs,(xfMs+L0)/L0,'1',c=col)
L0s=H0s/ARs
inds1=np.where(ARs>3)

fit2 = np.polynomial.polynomial.polyfit(np.log(ARs[inds1]), np.log((xfMs[inds1]+L0s[inds1])/(L0s[inds1])), [0, 1])

[line1]=ax3.plot(ARs_fit,np.exp(fit2[0])*ARs_fit**fit2[1], '-k',
        label=r"fit narrow: "+ format(np.exp(fit2[0]),'.1f') +r" $a^{"+ format(fit2[1],'.2f') +r"}$",lw=1)

ax3.legend()

#%%




# for sR in R_2d:
#     sR.scLength(0.01)
#     nF = int(sR.dConfig['nFrames'])
#     sR.nYplot=int(sR.resY//2)
#     sR.loadVTK(int(nF * sR.dConfig['fps'] / 15))


#%%
plt.close('all')
fig, ax = plt.subplots(1,1)
mus=[]
hs=[]
for sR in R_2d:
    hmax=[]
    for k in range(3,6):
        hmax.append(np.mean(sR.grid_y[0,np.where(sR.reshaped_Phi[k,:]>0.5)[0][-1]]))
    sR.hmax=np.mean(hmax)
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
    # mus=np.array([0.44, 0.45, 0.46, 0.48, 0.53, 0.58])
    # hs=np.array([0.05885417, 0.061679173, 0.0635625, 0.0678, 0.079100005, 0.08851667])
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
Ws=[1,10,6,15,2,20]

mu_eqs=[]
hmaxs=[]
i=0
for sR,w in zip(R_3d,Ws):
    H0=sR.dConfig['box'][2]*0.8
    hmax=[]
    for k in range(2,5):
        hmax.append(sR.grid_z[0,0,np.where(sR.reshaped_Phi[k,3,:]>0.5)[0][-1]])
    sR.hmax=np.mean(hmax)
    mu_eq=xopt.x[0]+xopt.x[1]*sR.hmax**xopt.x[2]
    # print(sR.hmax)
    mu_eqs.append(mu_eq)
    # ax.plot(w,sR.hmax,'k+')
    #evaluate the lost
    V = area(sR.findContourPhi(level=0.5)[0])
    hmaxs.append(sR.hmax)
    lost = (Vini[i]-V)/Vini[i]*100
    print(lost)
    i+=1

Ws=np.array(Ws)

mu_eq_err=xopt.x[1]*(sR.hmax+0.004)**xopt.x[2]-(xopt.x[1]*(sR.hmax-0.004)**xopt.x[2])

ax.plot((Ws*0.01/H0)**(-1),mu_eqs,'p', color='darkblue',label=r'3d simulation  $\mu=0.44$',ms=4)

ax.errorbar((Ws*0.01/H0)**(-1),mu_eqs, yerr=mu_eq_err, color='k',fmt='none', markersize=8,markeredgewidth=1, capsize=2, )

ax.set_ylabel(r'$\mu_{2D,eq}$')
ax.set_xlabel(r'$H_0/W$ ')

ws_Ls=np.linspace(0.01,0.30,100)
muIonescu = 0.38 + 0.18 * 0.05 / ws_Ls
ax.plot((ws_Ls/H0)**(-1),muIonescu,'--',c='darkred',label=' Ionescu et al. (2015) Correction',lw=1.4)



ax.set_position([0.15,0.15,0.8,0.8])
ax.legend()
plt.show()

fig.savefig(driveFolder+"/programs/gitLab/dry-granular/doc/article/figures/AR2_mu2Deq_W.pdf", dpi=150)





# Save data in csv
#%%
import csv 
filename = "/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/walls/AR2_0deg/AR2_0deg_3D.csv"
slope= np.absolute(sR.slope)
variables = [['Ws',Ws],['H0',np.ones(len(mu_eqs))*H0],['slope',np.ones(len(mu_eqs))*slope],['Hf',hmaxs],['mu_2D_eq',mu_eqs]]
opyf.write_csvScalar(filename, variables)


#%%
import opyf
import matplotlib.pyplot as plt
filename = "/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/walls/AR2_0deg/AR2_0deg_3D.csv"

header, data= opyf.read_csv(filename)

plt.rcParams["text.usetex"]=True
plt.close('all')
fig, ax = plt.subplots(1,1,figsize=(4,3))
Ws = data[:,0]
mu_eqs = data[:,4]
H0 = data[0,1]
ax.plot((Ws*0.01/H0)**(-1),mu_eqs,'p', color='darkblue',label=r'3d simulation  $\mu=0.44$',ms=4)

ax.errorbar((Ws*0.01/H0)**(-1),mu_eqs, yerr=0.02, color='k',fmt='none', markersize=8,markeredgewidth=1, capsize=2, )
ax.set_position([0.15,0.15,0.8,0.8])
ax.legend()
plt.show()
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
