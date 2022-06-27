
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% Author : Gauthier Rousseau
    # intialize exteral packages
import sys
import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
# sys.path.append('/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python')
import d6py
import matplotlib.pyplot as plt
# sys.path.append('/media/gauthier/Data-Gauthier/programs/gitHub/opyflow')
import opyf #from opyflow library some rendering function may be employed

    
# intialize font type and size
plt.rcParams['font.size'] = 8.0
plt.rcParams['font.family']= 'serif'
plt.rcParams['xtick.labelsize'] = 6.0
plt.rcParams['ytick.labelsize']=6.0
# print(r'\includegraphics{test_2.pdf}')
driveFolder='/media/gauthier/Data-Gauthier/Gauthier'
maind6OutFolder='/media/gauthier/Samsung_T5/sand6_out/'
paths,folders,listDictConf,listNumRun=d6py.findOutSand6Paths(maind6OutFolder,4)

N=7
# Reference
selectedRuns, selectedDict = d6py.whereSand6OutFromParms(listNumRun,mute=True,substeps=20,muRigid=0.18,mu=0.38, delta_mu=0, runNumber=N, dimSim=3)

#%%
scale = 0.01  #1cm


for sR in selectedRuns:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])
    # d6py.d62vtk(sR.d6OutFolder+'/')
    # sR.scLength(scale)

mainExpFolder=driveFolder+'/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs_opyf'
runExp1=d6py.ExperimentalRun(N,mainExpFolder,loadField=True)
runExp1.scLength(scale)

vec_cap_time=np.array([0,1,2,4,8,16,32,64,128])*runExp1.typicalTime

plt.ion()
plt.close('all')

# %matplotlib qt5
colors = [(33./255,66./255,99./255,0.1),(1,1,0.3,0.9), (0.8,0,0,0.9), (0,0,0,0.9)]
position = [0,0.1,0.5, 1]
cmap = opyf.make_cmap(colors, position)

#%%

cmap.set_over(color='g')
cmap.set_under(alpha=0)
L = runExp1.dictE['L'] / runExp1.scaleLength
H = (runExp1.dictE['H'] + 0.01) / runExp1.scaleLength

w_fig = (L + runExp1.xmax * runExp1.dictE['H'] / runExp1.scaleLength) / 83.4 * 8

plt.rcParams['font.size'] = 9



for ifile in [4]:
    # fig, [ax, ax1] = plt.subplots(2, 1)
    fig = plt.figure(dpi=142,figsize=(w_fig*0.39 ,3*0.39))
    # fig.set_size_inches(, forward=True)
    ax = plt.gca()
    Y_lim =[-0.01/runExp1.scaleLength,H]
    X_lim=[-L,runExp1.xmax*H]


    kt=0
    ax.cla()
    # ax1.cla()
    h1, l1 = [], []
    ls=['--','-','-.','--']
    c = [1.5, 1.3, 1.1, 0.9]
    NsR=len(selectedRuns)
    for sR,i in zip(selectedRuns,range(len(selectedRuns))):

        sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
        # sR.plotVelocity(ax,vmin=0,vmax=.5,cmap=cmap,interpolation='gaussian',alpha=0.5)
        # sR.plotI(ax,vmin=0,vmax=.1)
#        sR.plotPressure(ax,vmin=0,vmax=1000)
        # ax.imshow(np.mean(sR.P[:,:,:],axis=1),extent=sR.extentR,vmin=0,vmax=2000)
        # ax.imshow(np.mean(sR.epsilon21[:,:,:],axis=1),extent=sR.extentR,vmin=-50,vmax=50)
        sR.plotContour(ax,levels=[0.5],linewidths=c[i%4],linestyles=ls[i%4],colors=[cmap((NsR-i)/(len(selectedRuns)+2))]) 
        if i==0:
            # sR.plotI(ax,vmin=0,vmax=.1)
            # sR.opyfPointCloudScatter(ax)
            # im= sR.plotPhi(ax,cmap=cmap,vmin=0.9,vmax=2,interpolation='gaussian')
            # sR.plotPoints(ax,ls='',marker='.',markersize=0.2,color='k',alpha=0.7)
            im = sR.opyfPointCloudColoredScatter(ax,nvec=4000, mute=True, vmin=0, vmax=1, s=2, cmap=cmap, rasterized=True)
            # im_temp = im.copy()
            # im_temp[:,:, 0:3] = im[:,:, 1:4]
            # im_temp[:,:,3]= im[:,:, 0]
            # plt.ion()
            # plt.figure()
            # plt.imshow(im_temp)
            # plt.show()
            # sR.plotVelocity(ax,vmin=0,vmax=1,cmap=cmap,interpolation='gaussian')
        # if i==0:
            # sR.opyfPointCloudColoredScatter(ax,vmin=0,vmax=.5,s=4,cmap=cmap,alpha=0.5,nvec=3000)
            sR.plotDoor(ax, alpha=0.5)
            # sR.plotDoor(ax1,alpha=0.5)
            # except:
            #     print('no door')
 
        # h,l = sR.CS.legend_elements(str(sR.dimSim)+"D- \mu_{rb}="+format(sR.dConfig['muRigid'],'0.2f') +"D- frac_h="+format(sR.dConfig['frac_h'],'0.2f') + "- \mu= "+ format(sR.dConfig['mu'],'0.2f')+ "- delta_mu= "+ format(sR.dConfig['delta_mu'],'0.4f')+" - \phi")        
        # h,l = sR.CS.legend_elements(str(sR.dimSim)+"D-~\delta_t="+ toS(1000/sR.dConfig['substeps']/sR.dConfig['fps'],2)+ "~$ms$ ~-~ \mu= "+ toS(sR.dConfig['mu'],2)+ " - \phi")        
        # h,l = sR.CS.legend_elements(str(sR.dimSim)+"D-~\mu_{RB}="+ toS(sR.dConfig['muRigid'],2)+ " ~-~ \mu= "+ toS(sR.dConfig['mu'],2)+ " - \phi")        

        # h1, l1=h1+h, l1+l
        
    # def opyfColorBar(fig, im, label='Magnitude [px/Dt]', **args):

    # return fig, cb
    
    # runExp1.plotField(ax1,np.max([ifile*10-5,0]),vmin=0,vmax=1,cmap=cmap)
    # runExp1.plotDepthProfile(ax1,np.max([ifile*10-5,0]),linestyle='--',color='k',linewidth=1,label="Experience")
    # runExp1.plotDepthProfile(ax,np.max([ifile*10-5,0]),linestyle='--',color='k',linewidth=1,label="Experience")
    # ax.legend(h1+runExp1.l, l1+[r"Expe - tan($\theta$)=0.43 $\pm$ 0.03"],fontsize=10,framealpha=0,loc=1)
    # cbaxes = fig.add_axes([0.15, 0.1, 0.70, 0.03])

    # cb = fig.colorbar(im, cax=cbaxes, orientation='horizontal')
    # cb.set_label('Velocity [m/s]')
    
    ax.axis('equal')
    # ax1.axis('equal')
    ax.set_xlabel("X [cm]",fontsize=7)

    ax.set_ylabel("Y [cm]",fontsize=7)
    # ax1.set_ylabel("Y [cm]")
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    # ax1.set_xticklabels([])
    # ax.set_ylim([-2,12])
    # ax.set_xlim([-22.6,50])
    # ax1.set_ylim([-2,12])
    # ax1.set_xlim([-22.6,50])
    [x,y,X,Y]=ax.get_position().bounds
    ax.set_position([0.15, 0.02, 0.7, 0.7])
    # [x,y,X,Y]=ax1.get_position().bounds
    # ax1.set_position([0.1,y+0.05, X-0.1, Y])

    # d6py.setFigure(fig,ax1,runExp1.scaleLength,unit='cm')
#    d6py.setFigure(fig,ax1,runExp1.scaleLength,unit='cm')

    ax.plot([-L, -L], [-10, 2*H], '-k',linewidth=2)
    # ax1.plot([-L,-L],[-1,H],'-k',linewidth=6)
    # ax1.set_ylim(np.array([-0.01,runExp1.dictE['H']+0.01])/runExp1.scaleLength)
    # ax1.set_xlim(np.array([-runExp1.dictE['L'],runExp1.xmax*runExp1.dictE['H']])/runExp1.scaleLength)

    # ax1.text(28,6.5,s='Time='+format(runExp1.Time[np.max([ifile*10-5,0])],'10.2f')+' s', fontsize=8,zorder=2)

    # plt.show()
    # outFolder=listSaveFolders[runExp1.runNumber]+'/'+sR.dConfig['folder']
    # opyf.mkdir2(outFolder)
    # fig.savefig(outFolder+'/test_DeltaT_compMPM_EXP_t='+format(runExp1.Time[ifile],'02.2f')+'.svg')
    # plt.pause(0.1)
    ax.set_ylim(Y_lim)
    ax.set_xlim([X_lim[0], 30])
    plt.draw() 
    plt.show()
    plt.pause(0.2)

fig.savefig("test_savefig_pdf.pdf", dpi=300)
print(r'\includegraphics{test_savefig_pdf.pdf}')


# %%
