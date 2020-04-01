
#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Author : Gauthier Rousseau
    # intialize exteral packages
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
# sys.path.append('/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python')
import d6py
from d6py.Tools import *
import matplotlib.pyplot as plt
# sys.path.append('/media/gauthier/Data-Gauthier/programs/gitHub/opyflow')
import opyf #from opyflow library some rendering function may be employed

sys.stdout = open(os.devnull, 'w')    
# intialize font type and size
plt.rcParams['font.size'] = 8.0

plt.rcParams['xtick.labelsize'] = 8.0
plt.rcParams['ytick.labelsize'] = 8.0
plt.rcParams['ytick.labelsize'] = 8.0
plt.rcParams['axes.linewidth'] = 1
# print(r'\includegraphics{test_2.pdf}')
driveFolder='/media/gauthier/Data-Gauthier/Gauthier'
maind6OutFolder='/media/gauthier/Samsung_T5/sand6_out/'
paths,folders,listDictConf,listNumRun=d6py.findOutSand6Paths(maind6OutFolder,4)

N=7
# Reference
Ref, selectedDict = d6py.whereSand6OutFromParms(listNumRun,mute=True,substeps=20,muRigid=0.18,mu=0.38, delta_mu=0, runNumber=N, dimSim=3)
selectedRuns, selectedDict = d6py.whereSand6OutFromParms(listNumRun,mute=True,substeps=20,muRigid=0.18,mu=0.43, delta_mu=0, runNumber=N, dimSim=3)

selectedRuns=[Ref[0],selectedRuns[0]]


scale = 0.01  #1cm


for sR in selectedRuns:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])


mainExpFolder=driveFolder+'/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs_opyf'
runExp1=d6py.ExperimentalRun(N,mainExpFolder,loadField=True)
runExp1.scLength(scale)

runExp1.loadVideo('/media/gauthier/Samsung_T5/MPM_data/Collapse_experiment/Video_src',mute=True)
plt.close('all')
step,shift,Ntot=8,5,(int(runExp1.dictE['nFrames'])+20)*2
runExp1.video.set_vecTime(starting_frame=int(runExp1.dictE['framedeb']) - 100, step=step, shift=shift, Ntot=Ntot)
OR2 = runExp1.dictE['OR2']
sys.path.append('/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python/imageProcessing')
from Tools_collapses import mask_collapses

BW,vecXE,HE=mask_collapses(opyf.Tools.convertToGrayScale(runExp1.video.cropFrameInit),OR2,0.8,runExp1.dictE['scale'])
OR2[1]=OR2[1]-np.mean(HE[-int(len(HE)/2):-30])/runExp1.dictE['scale']


runExp1.video.scaleData(metersPerPx=runExp1.dictE['scale']/scale, framesPerSecond=runExp1.dictE['fps'], origin=OR2)
runExp1.video.paramPlot['extentFrame']
vec_cap_time=np.array([0,1,2,4,8,16,32,64,128])*runExp1.typicalTime
#%%
plt.ion()
plt.rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
plt.close('all')

cmapg= opyf.custom_cmap.make_cmap_customized(Palette='green')

colors = [(33./255,66./255,99./255,0.05),(1,1,0.3,0.9), (0.8,0,0,0.9), (0,0,0,0.9)]
position = [0,0.1,0.5, 1]
cmap = opyf.make_cmap(colors, position)



# cmap.set_over(color='g')
cmap.set_under(alpha=0)
L = runExp1.dictE['L'] / runExp1.scaleLength
H = (runExp1.dictE['H'] + 0.01) / runExp1.scaleLength

w_fig = (L + runExp1.xmax * runExp1.dictE['H'] / runExp1.scaleLength) / 83.4 * 8

w_fig= 15 
N = 3
M = 3
fig, axs = plt.subplots(N, M, dpi=142, figsize=(w_fig * 0.39, 13 * 0.39))

# fig = plt.figure(dpi=142, figsize=(w_fig * 0.39, 11 * 0.39))

Y_lim =[-0.01/runExp1.scaleLength,H]
X_lim = [-L-1, runExp1.xmax * H]
w_axs = 0.28
h_axs = 0.2
w_s = 0.09
w_s2 =0.02
for i in range(N):
    for j in range(M):

        axs[i, j].plot([-L, -L], [-10, 2 * H], '-k', linewidth=2)
        [x, y, X, Y] = axs[i, j].get_position().bounds
        axs[i, j].set_position([w_s+(w_axs + w_s2) * j, 0.2+(h_axs - 0.03) * (N-i), w_axs, h_axs])
        axs[i, j].set_xlim([X_lim[0], 7])
        axs[i, j].set_ylim(Y_lim)
        axs[i, j].set_aspect('equal',adjustable='box')
        axs[i, j].set_yticklabels([])
        axs[i, j].set_xticklabels([])

#draw lines to separate experimental plots
ax_draw = fig.add_axes([0, 0, 1, 1], zorder=-10)
ax_draw.grid()
X_line = w_s + (w_axs + 0.01)
[x, y, X, Y] = axs[2, 2].get_position().bounds 
Y_sep = y-0.03
ax_draw.plot([X_line, X_line], [Y_sep, 0.98], linewidth=0.5, color='k')
ax_draw.plot([0.02, 0.98], [Y_sep, Y_sep], linewidth=0.5, color='k')

# Add a big axe to draw the final state entirely
ax2 = fig.add_axes([w_s, 0.00, x + X -w_s, Y_sep-0.03], zorder=-10)
ax2.set_yticklabels([])
ax2.set_xticklabels([])
ax2.set_xlim([X_lim[0], 52])
ax2.plot([-L, -L], [-10, 2 * H], '-k', linewidth=2)
ax2.set_ylim(Y_lim)
ax2.set_aspect('equal',adjustable='box')
# draw scale and axes


for i in range(M):
    for j in range(M):
        axs[i, j].plot([-200, 200], [0, 0], 'k', linewidth=0.5)
        axs[i, j].plot([0, 0], [-200, 200], 'k', linewidth=0.5)
        axs[i, j].text(-0.5, 0, 'O', fontsize=8, horizontalalignment='right', verticalalignment='bottom', bbox=dict(boxstyle="round", fc="w", alpha=0.2))
        draw_sc=0
        if i == 2:
            draw_sc=1
            x_sc, y_sc = -15, 1
            lx_sc, ly_sc = 10, 3.8
            if j == 2: ly_sc=4.3     
        if i == 0 and j == 0:
            draw_sc=1
            x_sc, y_sc = -15, 1
            lx_sc, ly_sc = 10, 5
        if i == 0 and j == 1:
            x_sc, y_sc = -17.5, 2.5
            lx_sc, ly_sc = 6, 6
            arrowprops = dict(
            arrowstyle = "<-",
             color='black',
             linewidth=1)
            axs[i,j].annotate('', xy=(x_sc-1, y_sc), xytext=(x_sc+lx_sc, y_sc),
            arrowprops=arrowprops, fontsize=6, zorder=20)
            axs[i,j].annotate('', xy=(x_sc, y_sc-1), xytext=(x_sc, y_sc+ly_sc),
            arrowprops=arrowprops, fontsize=6, zorder=20)
            axs[i,j].text(x_sc+lx_sc, y_sc,'$x$',fontsize=8, verticalalignment='bottom', horizontalalignment='center')
            axs[i,j].text(x_sc, y_sc+ly_sc,'$y$', verticalalignment='center', horizontalalignment='right',fontsize=8)
        if draw_sc==1:    
            axs[i, j].plot([x_sc, x_sc+lx_sc], [y_sc, y_sc], 'k', linewidth=1)
            axs[i, j].plot([x_sc, x_sc], [y_sc, y_sc + ly_sc], 'k', linewidth=1)
            axs[i, j].text(x_sc + lx_sc / 2, y_sc + 1, format(lx_sc, '.1f') + ' cm', fontsize=6, horizontalalignment='center', bbox=dict(boxstyle="round", fc="w", alpha=0.5))
            axs[i, j].text(x_sc - 1, y_sc + ly_sc / 2, format(ly_sc, '.1f') + ' cm', fontsize=6, horizontalalignment='right', verticalalignment='center', bbox=dict(boxstyle="round", fc="w", alpha=0.5))
                   
ax2.plot([-200, 200], [0, 0], 'k', linewidth=0.5)
ax2.plot([0, 0], [-200, 200], 'k', linewidth=0.5)
ax2.text(-0.5, 0, 'O', fontsize=8, horizontalalignment='right', verticalalignment='bottom', bbox=dict(boxstyle="round", fc="w", alpha=0.5))

#draw gravity vector
x_sc, y_sc = -20, 9
lx_sc, ly_sc = 10, 5
vecgy=y_sc
vecgx=x_sc
lg=8
dxT=0.2*lg
dyT = 0.1 * lg
slopeinrad = 15 * 3.15 / 180
for j in [2]:
    axs[0,j].quiver(vecgx, vecgy,lg*np.sin(slopeinrad),-lg*np.cos(slopeinrad),width=0.008,linewidth=1,angles='xy', scale_units='xy', scale=1)
    axs[0,j].plot([vecgx, vecgx],[vecgy,vecgy-1.2*lg],'k',linewidth=1)
    from matplotlib.patches import Circle, Wedge, Polygon, Arc
    arc=Arc([vecgx,vecgy],lg*1.4,lg*1.4,theta1=270,theta2=270+np.rad2deg(slopeinrad),linewidth=0.5,color='k')
    #    arc=Arc([vecgx,vecgy],lg*1.4,lg*1.4,theta1=270,theta2=270+np.rad2deg(slopeinrad),linewidth=2,color='k')
    #  
    axs[0,j].add_patch(arc)
    axs[0,j].text(vecgx+lg*slopeinrad+dxT,vecgy-lg*np.cos(slopeinrad),r'$\vec{g}$',fontsize=7)
    axs[0,j].text(vecgx+lg*slopeinrad+dxT,vecgy-lg*np.cos(slopeinrad)+5*dyT,r'$i$='+format(15,'1.0f')+r'$^\circ $',fontsize=7)    



# write times
[x, y, X, Y] = axs[0, 0].get_position().bounds    
plt.figtext(0.01,y+Y/2,'t=0 s', fontsize=7)
[x, y, X, Y] = axs[1, 0].get_position().bounds    
plt.figtext(0.01, y + Y / 2, 't=0.2 s', fontsize=7)
[x, y, X, Y] = axs[2, 0].get_position().bounds  
plt.figtext(0.01, y + Y / 2, 't=0.6 s', fontsize=7)
[x, y, X, Y] = ax2.get_position().bounds  
plt.figtext(0.01,y+Y/2,'t=1.6 s', fontsize=7)






# write Experiments and Model
[x, y, X, Y] = axs[0, 0].get_position().bounds    
plt.figtext(x+X/2,y+Y+0.07,'Experiment', fontsize=10,horizontalalignment='center')
[x, y, X, Y] = axs[0, 1].get_position().bounds    
plt.figtext(x+X+w_s2/2,y+Y+0.07,'Simulations', fontsize=10,horizontalalignment='center')
plt.figtext(x + X / 2, y + Y + 0.03, '$\mu_1$=0.38', fontsize=10, horizontalalignment='center')
[x, y, X, Y] = axs[0, 2].get_position().bounds 
plt.figtext(x+X/2,y+Y+0.03,'$\mu_1$=0.43', fontsize=10,horizontalalignment='center')

ax_draw.set_xlim([0, 1])
ax_draw.set_ylim([0, 1])


arrowprops = dict(
    arrowstyle = "-",
    # connectionstyle="angle,angleA=0,angleB=90,rad=2",
    # shrink=0.05,
    # width=1.,
    # headwidth=4.,
    # headlength=1,
    color='black')
    

axs[1,0].annotate('uplifting gate', xy=(2, 15), xytext=(-3, 14.5),
            arrowprops=arrowprops, fontsize=6, zorder=20, verticalalignment='center', horizontalalignment='right',bbox=dict(boxstyle="round", fc="w",alpha=0.5))

arrowprops = dict(
    arrowstyle = "->",
    color='black')

axs[1,0].annotate('', xy=(4, 15), xytext=(4, 10),
            arrowprops=arrowprops, fontsize=6, zorder=20, verticalalignment='center', horizontalalignment='right')            
           
# for i in range(N):

#         [x, y, X, Y] = axs[i, 1].get_position().bounds
#         axs[i, 1].set_position([0.31, y, X, Y])
# axs[1,0].set_ylim([-20,20])          
plt.show()

# %
k=0
for ifile in [0, 3, 9]:
    
    ax=axs[k,0]
    h1, l1 = [], []
    ls=['--','-','-.','--']
    c = [1., 1.2, 1.1, 0.9]
    NsR = len(selectedRuns)
    
    runExp1.video.readFrame(int(runExp1.dictE['framedeb']) - 100 + 10 * 240)
    

    
    runExp1.video.readFrame(np.max([int(runExp1.dictE['framedeb']) - 100 + ifile * (100) - 50, int(runExp1.dictE['framedeb']) - 100]))
    vis=opyf.Render.CLAHEbrightness(runExp1.video.vis,140)
    ax.imshow(vis, extent=runExp1.video.paramPlot['extentFrame'])
    runExp1.plotField(ax, np.max([ifile * 10 - 5, 0]), vmin=0, vmax=1, cmap=cmap)
    runExp1.plotDepthProfile(ax, np.max([ifile * 10 - 5, 0]), linestyle='-.', color='k', linewidth=1, label="Experience")
    # axs[k,0].set_ylabel("Y [cm]",fontsize=7)

    for sR, i in zip(selectedRuns, range(len(selectedRuns))):
        ax = axs[k, i + 1]
        sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
        # sR.plotVelocity(ax,vmin=0,vmax=.5,cmap=cmap,interpolation='gaussian',alpha=0.5)
        # sR.plotI(ax,vmin=0,vmax=.1)
#        sR.plotPressure(ax,vmin=0,vmax=1000)
        # ax.imshow(np.mean(sR.P[:,:,:],axis=1),extent=sR.extentR,vmin=0,vmax=2000)
        # ax.imshow(np.mean(sR.epsilon21[:,:,:],axis=1),extent=sR.extentR,vmin=-50,vmax=50)

        sR.plotContour(ax,levels=[0.5],linewidths=c[i%4],linestyles=ls[i%4],colors=[cmapg((NsR-i)/(NsR+2))]) 
            # sR.plotI(ax,vmin=0,vmax=.1)
            # im= sR.plotPhi(ax,cmap=cmap,vmin=0.9,vmax=2,interpolation='gaussian')
            # sR.plotPoints(ax,ls='',marker='.',markersize=0.2,color='k',alpha=0.7)
        im = sR.opyfPointCloudColoredScatter(ax,nvec=5000, mute=True, vmin=0, vmax=1, s=0.8, cmap=cmap, rasterized=True)
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
        runExp1.plotDepthProfile(ax,np.max([ifile*10-5,0]),linestyle='-.',color='k',linewidth=1,label="Experience")
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        if ifile == nF:
            sR.plotContour(ax2,levels=[0.5],linewidths=c[i%4],linestyles=ls[i%4],colors=[cmapg((NsR-i)/(NsR+2))]) 
            
            # sR.plotDoor(ax1,alpha=0.5)
            # except:
            #     print('no door')
 
        # h,l = sR.CS.legend_elements(str(sR.dimSim)+"D- \mu_{rb}="+format(sR.dConfig['muRigid'],'0.2f') +"D- frac_h="+format(sR.dConfig['frac_h'],'0.2f') + "- \mu= "+ format(sR.dConfig['mu'],'0.2f')+ "- delta_mu= "+ format(sR.dConfig['delta_mu'],'0.4f')+" - \phi")        
        # h,l = sR.CS.legend_elements(str(sR.dimSim)+"D-~\delta_t="+ toS(1000/sR.dConfig['substeps']/sR.dConfig['fps'],2)+ "~$ms$ ~-~ \mu= "+ toS(sR.dConfig['mu'],2)+ " - \phi")        
        h= sR.CS.legend_elements(str(sR.dimSim)+"D-~\mu_{RB}="+ toS(sR.dConfig['muRigid'],2)+ " ~-~ \mu= "+ toS(sR.dConfig['mu'],2)+ " - \phi")[0]        
        l=[str(sR.dimSim)+r"$D-~\mu_{RB}$="+ toS(sR.dConfig['muRigid'],2)+ r"$ ~-~ \mu$= "+ toS(sR.dConfig['mu'],2)]
        h1, l1=h1+h, l1+l
        


        # runExp1.plotDepthProfile(ax1,np.max([ifile*10-5,0]),linestyle='--',color='k',linewidth=1,label="Experience")
        # ax.legend(h1+runExp1.l, l1+[r"Expe - tan($\theta$)=0.43 $\pm$ 0.03"],fontsize=10,framealpha=0,loc=1)

        
        # ax1.set_xticklabels([])
        # ax.set_ylim([-2,12])
        # ax.set_xlim([-22.6,50])
        # ax1.set_ylim([-2,12])
        # ax1.set_xlim([-22.6,50])
        # [x,y,X,Y]=ax.get_position().bounds
        # ax.set_position([0.1, 0.02, 0.8, 0.8])
        # [x,y,X,Y]=ax1.get_position().bounds
        # ax1.set_position([0.1,y+0.05, X-0.1, Y])

        # d6py.setFigure(fig,ax1,runExp1.scaleLength,unit='cm')
    #    d6py.setFigure(fig,ax1,runExp1.scaleLength,unit='cm')


        # ax1.plot([-L,-L],[-1,H],'-k',linewidth=6)
        # ax1.set_ylim(np.array([-0.01,runExp1.dictE['H']+0.01])/runExp1.scaleLength)
        # ax1.set_xlim(np.array([-runExp1.dictE['L'],runExp1.xmax*runExp1.dictE['H']])/runExp1.scaleLength)

        # ax1.text(28,6.5,s='Time='+format(runExp1.Time[np.max([ifile*10-5,0])],'10.2f')+' s', fontsize=8,zorder=2)

        # plt.show()
        # outFolder=listSaveFolders[runExp1.runNumber]+'/'+sR.dConfig['folder']
        # opyf.mkdir2(outFolder)
        # fig.savefig(outFolder+'/test_DeltaT_compMPM_EXP_t='+format(runExp1.Time[ifile],'02.2f')+'.svg')

        # plt.draw() 

    ax = axs[k, 0]
    sR.plotDoor(ax, alpha=0.5)
        # plt.pause(0.2)
    k += 1



ifile = nF
h1, l1 = [], []
runExp1.video.readFrame(np.max([int(runExp1.dictE['framedeb']) - 100 + ifile * (100) + 50, int(runExp1.dictE['framedeb']) - 100]))
vis=opyf.Render.CLAHEbrightness(runExp1.video.vis,40)
ax2.imshow(vis, extent=runExp1.video.paramPlot['extentFrame'])
# runExp1.plotField(ax2, np.max([ifile * 10 - 5, 0]), vmin=0, vmax=1, cmap=cmap)
runExp1.plotDepthProfile(ax2, np.max([ifile * 10 - 5, 0]), linestyle='-.', color='k', linewidth=1, label="Experience")
for sR, i in zip(selectedRuns, range(len(selectedRuns))):
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    sR.plotContour(ax2,levels=[0.5],linewidths=c[i%4],linestyles=ls[i%4],colors=[cmapg((NsR-i)/(NsR+2))])         
    h= sR.CS.legend_elements(str(sR.dimSim)+"D-~\mu_{RB}="+ toS(sR.dConfig['muRigid'],2)+ " ~-~ \mu= "+ toS(sR.dConfig['mu'],2)+ " - \phi")[0]        
    l=[str(sR.dimSim)+ r"D-~$ ~ \mu_1$= "+ toS(sR.dConfig['mu'],2)+r"$~-~\delta_{\mu}$="+ toS(sR.dConfig['delta_mu'],2)]
    h1, l1=h1+h, l1+l
ax2.legend(h1 + runExp1.l, l1 + [r"Exp. - tan($\theta$)=0.43 $\pm$ 0.03"], fontsize=7, framealpha=0.5, loc=1)


[x, y, X, Y] = axs[2, 2].get_position().bounds 
cbaxes = fig.add_axes([0.18, y-0.08, 0.70, 0.03])

cb = fig.colorbar(im, cax=cbaxes, orientation='horizontal',extend='both')
cb.set_label('Velocity norm [m/s]', fontsize=8)


[x, y, X, Y] = ax2.get_position().bounds 
cbaxes = ax2.set_position([x, y-0.03, X, Y])
plt.show()
# fig.savefig("test_savefig_pdf.pdf", dpi=300)
# fig.savefig("test_savefig_pdf.svg", dpi=150)
sys.stdout = sys.__stdout__
print(r'\includegraphics{test_savefig_pdf.pdf}')


# %%
