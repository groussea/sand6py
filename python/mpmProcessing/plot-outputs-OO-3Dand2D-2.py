#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:52:28 2019

@author: Gauthier Rousseau
"""


#%%

#
# 

# %matplotlib inline

import sys 
sys.path.append('/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python')
import d6py
from d6py.Tools import *
sys.path.append('/media/gauthier/Data-Gauthier/programs/gitHub/opyflow')
import opyf
sys.path

import numpy as np
import matplotlib.pyplot as plt
# driveFolder='/scratch/garousse/'
# driveFolder='/media/gauthier/Gauthier_Backup/'
# driveFolder='/media/garousse/Gauthier_Backup/'
driveFolder='/media/gauthier/Data-Gauthier/Gauthier'
maind6OutFolder=driveFolder+'/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs/Tests/'
maind6OutFolder='/media/gauthier/Samsung_T5/sand6_out/'
paths,folders,listDictConf,listNumRun=d6py.findOutSand6Paths(maind6OutFolder,4)



#%

N=7

selectedRuns, selectedDict = d6py.whereSand6OutFromParms(listNumRun,muRigid=0.18,mu=0.38, delta_mu=0, runNumber=N, dimSim=3)





selectedRuns = [selectedRuns[2],   selectedRuns[0],selectedRuns[3],  selectedRuns[1]]


# selectedRuns = [selectedRuns[1], selectedRuns[2]] # 2D
#%%
for sR in selectedRuns:

    sR.scLength(0.01)
    print(sR.dConfig['substeps'])
    print(sR.dConfig['res'])
    print(sR.dConfig['mu'])
    nF = int(sR.dConfig['nFrames'])
    # d6py.d62vtk(sR.d6OutFolder+'/')
mainExpFolder=driveFolder+'/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs_opyf'
runExp1=d6py.ExperimentalRun(N,mainExpFolder,loadField=True)

runExp1.scLength(0.01)

plt.rc('text', usetex=True)#from opyflow library some rendering function may be employed

vec_cap_time=np.array([0,1,2,4,8,16,32,64,128])*runExp1.typicalTime

mainSaveFolder='/media/gauthier/Samsung_T5/plots_d6'
listSaveFolders=np.sort([mainSaveFolder+d for d in runExp1.dictExp])
for l in listSaveFolders:
    opyf.mkdir2(l)
cmap=opyf.make_cmap_customized_asym(Palette='asym_mountain_full')
#%
plt.close('all')
plt.ion()
%matplotlib qt5
colors = [(33./255,66./255,99./255,0.1),(1,1,0.3,0.9), (0.8,0,0,0.9), (0,0,0,0.9)]
position = [0,0.1,0.5, 1]
cmap = opyf.make_cmap(colors, position)


# cmap=plt.get_cmap('jet')
#d6Out_title= ''.join([(s!='_')*s for s in d6Out])
#ax.set_title(d6Out_title,fontsize=9)
# cmap=opyf.make_cmap_customized_asym(Palette='asym_mountain_full')
# for ifile in range(1, 22):
cmap.set_over(color='g')
cmap.set_under(alpha=0)
L = runExp1.dictE['L'] / runExp1.scaleLength
H = (runExp1.dictE['H'] + 0.01) / runExp1.scaleLength

w_fig = (L + runExp1.xmax * runExp1.dictE['H'] / runExp1.scaleLength) / 60 * 8

for ifile in [0,3,nF]:
    fig, [ax, ax1] = plt.subplots(2, 1)
    fig.set_size_inches((w_fig ,4.5))
    kt=0
    ax.cla()
    ax1.cla()
    h1, l1 = [], []
    ls=['--','-','-.','--']
    c = [1.8, 1.2, 1.4, 1.5]
    NsR=len(selectedRuns)
    for sR,i in zip(selectedRuns,range(len(selectedRuns))):

        sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))

        # sR.generateDepthProfile()
        # sR.plotVelocity(ax,vmin=0,vmax=.5,cmap=cmap,interpolation='gaussian',alpha=0.5)
        # sR.plotI(ax,vmin=0,vmax=.1)
#        sR.plotPressure(ax,vmin=0,vmax=1000)
        # ax.imshow(np.mean(sR.P[:,:,:],axis=1),extent=sR.extentR,vmin=0,vmax=2000)
        # ax.imshow(np.mean(sR.epsilon21[:,:,:],axis=1),extent=sR.extentR,vmin=-50,vmax=50)
        sR.plotContour(ax,levels=[0.5],linewidths=c[i%4],linestyles=ls[i%4],colors=[cmap((NsR-i)/(len(selectedRuns)+2))]) 
        if i==0:

            # sR.plotI(ax,vmin=0,vmax=.1)
            # sR.opyfPointCloudScatter(ax)
            #
            # im= sR.plotPhi(ax,cmap=cmap,vmin=0.9,vmax=2,interpolation='gaussian')
            # sR.plotPoints(ax,ls='',marker='.',markersize=0.2,color='k',alpha=0.7)
            sR.opyfPointCloudColoredScatter(ax,vmin=0,vmax=0.5,s=4,cmap=cmap,nvec=5000)
            # sR.plotVelocity(ax,vmin=0,vmax=1,cmap=cmap,interpolation='gaussian')
        # if i==0:
            # sR.opyfPointCloudColoredScatter(ax,vmin=0,vmax=.5,s=4,cmap=cmap,alpha=0.5,nvec=3000)
            sR.plotDoor(ax, alpha=0.5)
            sR.plotDoor(ax1,alpha=0.5)
            # except:
            #     print('no door')
            
        # h,l = sR.CS.legend_elements(str(sR.dimSim)+"D- \mu_{rb}="+format(sR.dConfig['muRigid'],'0.2f') +"D- frac_h="+format(sR.dConfig['frac_h'],'0.2f') + "- \mu= "+ format(sR.dConfig['mu'],'0.2f')+ "- delta_mu= "+ format(sR.dConfig['delta_mu'],'0.4f')+" - \phi")        
        h,l = sR.CS.legend_elements(str(sR.dimSim)+"D-~\delta_t="+ toS(1000/sR.dConfig['substeps']/sR.dConfig['fps'],2)+ "~$ms$ ~-~ \mu= "+ toS(sR.dConfig['mu'],2)+ " - \phi")        
        # h,l = sR.CS.legend_elements(str(sR.dimSim)+"D-~\mu_{RB}="+ toS(sR.dConfig['muRigid'],2)+ " ~-~ \mu= "+ toS(sR.dConfig['mu'],2)+ " - \phi")        

        h1, l1=h1+h, l1+l
        
    # def opyfColorBar(fig, im, label='Magnitude [px/Dt]', **args):

    # return fig, cb
    
    runExp1.plotField(ax1,np.max([ifile*10-5,0]),vmin=0,vmax=1,cmap=cmap)
    runExp1.plotDepthProfile(ax1,np.max([ifile*10-5,0]),linestyle='--',color='k',linewidth=1,label="Experience")
    runExp1.plotDepthProfile(ax,np.max([ifile*10-5,0]),linestyle='--',color='k',linewidth=1,label="Experience")
    ax.legend(h1+runExp1.l, l1+[r"Expe - tan($\theta$)=0.43 $\pm$ 0.03"],fontsize=10,framealpha=0,loc=1)
    cbaxes = fig.add_axes([0.15, 0.1, 0.70, 0.03])

    cb = fig.colorbar(runExp1.im, cax=cbaxes, orientation='horizontal')
    cb.set_label('Velocity [m/s]')
    
    ax.axis('equal')
    ax1.axis('equal')
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Y [cm]")
    ax1.set_ylabel("Y [cm]")
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax1.set_xticklabels([])
    # ax.set_ylim([-2,12])
    # ax.set_xlim([-22.6,50])
    # ax1.set_ylim([-2,12])
    # ax1.set_xlim([-22.6,50])
    [x,y,X,Y]=ax.get_position().bounds
    ax.set_position([0.1, y+0.02, X-0.1, Y])
    [x,y,X,Y]=ax1.get_position().bounds
    ax1.set_position([0.1,y+0.05, X-0.1, Y])

    # d6py.setFigure(fig,ax1,runExp1.scaleLength,unit='cm')
#    d6py.setFigure(fig,ax1,runExp1.scaleLength,unit='cm')

    ax.set_ylim(np.array([-0.01,runExp1.dictE['H']+0.01])/runExp1.scaleLength)
    ax.set_xlim(np.array([-runExp1.dictE['L'],runExp1.xmax*runExp1.dictE['H']])/runExp1.scaleLength)
    ax.plot([-L, -L], [-1, H], '-k',linewidth=6)
    ax1.plot([-L,-L],[-1,H],'-k',linewidth=6)
    ax1.set_ylim(np.array([-0.01,runExp1.dictE['H']+0.01])/runExp1.scaleLength)
    ax1.set_xlim(np.array([-runExp1.dictE['L'],runExp1.xmax*runExp1.dictE['H']])/runExp1.scaleLength)

    ax1.text(28,6.5,s='Time='+format(runExp1.Time[np.max([ifile*10-5,0])],'10.2f')+' s', fontsize=8,zorder=2)

    plt.show()
    outFolder=listSaveFolders[runExp1.runNumber]+'/'+sR.dConfig['folder']
    opyf.mkdir2(outFolder)
    fig.savefig(outFolder+'/test_DeltaT_compMPM_EXP_t='+format(runExp1.Time[ifile],'02.2f')+'.svg')
    # plt.pause(0.5)
#%%    
#
A=selectedRuns1[0]





from pyevtk.hl import pointsToVTK

path='file.vtu'

# Example 2

x=np.ascontiguousarray(A.pointsp[:,0],dtype=np.float32)
y=np.ascontiguousarray(A.pointsp[:,1],dtype=np.float32)
z=np.ascontiguousarray(A.pointsp[:,2],dtype=np.float32)
vol=np.ascontiguousarray(A.vol,dtype=np.float32)


pointsToVTK("./line_points", x, y, z, data = {"vol" : vol*1e8})
# assert (x.size == y.size == z.size)

#%%


npoints = x.size

# create some temporary arrays to write grid topology
offsets = np.arange(start = 1, stop = npoints + 1, dtype = 'int32')   # index of last node in each cell
connectivity = np.arange(npoints, dtype = 'int32')                    # each point is only connected to itself
cell_types = np.empty(npoints, dtype = 'uint8') 
   
cell_types[:] = vtk.vtkVertex.tid

w = VtkFile(path, VtkUnstructuredGrid)
w.openGrid()
w.openPiece(ncells = npoints, npoints = npoints)

w.openElement("Points")
w.addData("points", (x,y,z))
w.closeElement("Points")
w.openElement("Cells")
w.addData("connectivity", connectivity)
w.addData("offsets", offsets)
w.addData("types", cell_types)
w.closeElement("Cells")

_addDataToFile(w, cellData = None, pointData = data)

w.closePiece()
w.closeGrid()
w.appendData( (x,y,z) )
w.appendData(connectivity).appendData(offsets).appendData(cell_types)

_appendDataToFile(w, cellData = None, pointData = data)

w.save()



#%%

#    if runExp1.Time[ifile]>=vec_cap_time[kt]:
#        kt+=1
    outFolder=listSaveFolders[runExp1.runNumber]+'/'+sR.dConfig['folder']
    outFolder='/media/gauthier/Gauthier_Backup/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular/doc/article/images/texplots/data/Run_07'
    opyf.mkdir2(outFolder)
    folderExp=outFolder+'/experience'
    opyf.mkdir2(folderExp)
    variables=[['X',runExp1.vecXexpD],['h_surf',runExp1.expD[np.max(ifile*10-5,0)]]]
    opyf.write_csvScalar(folderExp+'/outExp-'+str(ifile)+'.csv',variables )

    folderSand6=outFolder+'/sand6'
    opyf.mkdir2(folderSand6)
    
    folderSand6=outFolder+'/sand6'
    opyf.mkdir2(folderSand6)   
 
    folderSand61=folderSand6+'/muRigid=0'
    opyf.mkdir2(folderSand61)  
    
    variables=[['X',selectedRuns[0].grid_x[:,0,0]],['h_surf',selectedRuns[0].h]]
    opyf.write_csvScalar(folderSand61+'/outSand6-'+str(ifile)+'.csv',variables )
    
    
    folderSand62=folderSand6+'/reference'
    opyf.mkdir2(folderSand62)   
    variables=[['X',selectedRuns[1].grid_x[:,0,0]],['h_surf',selectedRuns[1].h]]
    opyf.write_csvScalar(folderSand62+'/outSand6-'+str(ifile)+'.csv',variables )

    folderSand62=folderSand6+'/muRigid=0.5'
    opyf.mkdir2(folderSand62)   
    variables=[['X',selectedRuns[1].grid_x[:,0,0]],['h_surf',selectedRuns[1].h]]
    opyf.write_csvScalar(folderSand62+'/outSand6-'+str(ifile)+'.csv',variables )
        

    
    # instead of saving the datas we want to generate csv files for plottings with pgfplot
    
    
    # fig.savefig(outFolder+'/test2_compMPM_EXP_t='+format(runExp1.Time[ifile],'02.2f')+'.png')
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
A=sR.pointsp
sR.pointsLayer=[]
Ymin =  sR.grid_y[0,sR.nYplot,0]-sR.dy/2
            Ymax =  sR.grid_y[0,sR.nYplot,0]+sR.dy/2
for p in sR.pointsp:
    if p[1]> Ymin and p[1] < Ymax:
        sR.pointsLayer.append(p)
sR.pointsLayer=np.array(sR.pointsLayer)
ax.plot(sR.pointsLayer[:,0]/sR.scaleLength,sR.pointsLayer[:,2]/sR.scaleLength,**args)
