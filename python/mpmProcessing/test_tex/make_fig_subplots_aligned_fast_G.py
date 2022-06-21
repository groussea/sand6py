
# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Author : Gauthier Rousseau
import opyf  # from opyflow library some rendering function may be employed

sys.path.append(
    './..')
sys.path.append(
    './../..')
import d6py
import sys
sys.path.append(
    './../../imageProcessing')

from d6py.Tools import *
import sys
import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")
# plt.ioff()

paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(
    maind6OutFolder, 4)

runExps=[]
selectedRuns=[]
runExps=[]
selectedRuns=[]

for Nrun in range(4): # init
    scale = 0.01  # 1cm
    mainExpFolder = driveFolder + \
        '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'

    runExps.append(d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True))
    runExps[-1].scLength(scale)
    if Nrun < 4:
        mu=0.75
    else:
        mu=0.44
        
    R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, mu=mu, delta_mu=0., runNumber=Nrun, dimSim=3, delta_mu_start=0, keyWord='try2')

    selectedRuns.append(R1[0])
    
for sR in selectedRuns:
    sR.scLength(0.01)
    

#%%
import importlib
importlib.reload(plt)
#%%
plt.close('all')
plt.rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
cmapg = opyf.make_cmap_customized(Palette='green')

colors = [(33./255, 66./255, 99./255, 0.05),
        (1, 1, 0.3, 0.9), (0.8, 0, 0, 0.9), (0, 0, 0, 0.9)]
position = [0, 0.1, 0.5, 1]
cmap = opyf.make_cmap(colors, position)

cmap.set_under(alpha=0)


N = 4  # lignes
M = 3 # colonnes
fig, axs = plt.subplots(N, M, dpi=142, figsize=(7 , 4))

w_axs = 0.23
h_axs = 0.18
w_s = 0.02
w_s2 = 0.01
Largeur=30
Largeurf=75
for i in range(N):
    L = runExps[i].dictE['L'] / runExps[i].scaleLength-1
    H = (runExps[i].dictE['H']+0.02 ) / runExps[i].scaleLength
    L = runExps[0].dictE['L'] / runExps[0].scaleLength-1
    H=15
    Y_lim = [-0.015/runExps[i].scaleLength, H]
    X_lim = [-L-0.2, runExps[i].xmax * H] 
    for j in range(M):
        axs[i, j].plot([-L, -L], [-10, 2 * H], '-k', linewidth=2)
        axs[i, j].set_anchor('SW')
        [x, y, X, Y] = axs[i, j].get_position().bounds
        axs[i, j].set_position(
            [w_s+(w_axs + w_s2) * j, (h_axs+0.05 ) * (N-i-1)+0.06, w_axs*4, h_axs])
        axs[i, j].set_xlim([X_lim[0], Largeur+X_lim[0]])
        axs[i, j].set_ylim(Y_lim)
        axs[i, j].set_aspect('equal')
        axs[i, j].set_yticklabels([])
        axs[i, j].set_xticklabels([])
    axs[i,-1].set_xlim([X_lim[0], Largeurf+X_lim[0]])
# draw lines to separate experimental plots
ax_draw = fig.add_axes([0, 0, 1, 1], zorder=-10)
ax_draw.grid()
ax_draw.set_axis_off()
X_line = w_s + (w_axs + 0.01)

# draw scale and axes

my_bbox = dict(fc="w", alpha=0.3)

def draw_ax_sc(ax,x_sc, y_sc,lx_sc, ly_sc,fmt='.0f',shift_y_txt=0.5):
    ax.plot([x_sc, x_sc+lx_sc],
                    [y_sc, y_sc], '-+k', linewidth=1)
    ax.plot([x_sc, x_sc], [y_sc, y_sc + ly_sc],
                    '-+k', linewidth=1)
    ax.text(x_sc + lx_sc / 2, y_sc + shift_y_txt, format(lx_sc, fmt) + ' cm', fontsize=6,
                    horizontalalignment='center', bbox=my_bbox)
    # ax.text(x_sc - 1, y_sc + ly_sc / 2, format(ly_sc, fmt) + ' cm', fontsize=6,
    #                 horizontalalignment='right', verticalalignment='center', bbox=my_bbox)
    lx_sc, ly_sc = lx_sc +3, ly_sc+3
    arrowprops = dict(
        arrowstyle="<-",
        color='black',
        linewidth=1)
    ax.annotate('', xy=(x_sc-1, y_sc), xytext=(x_sc+lx_sc, y_sc),
                        arrowprops=arrowprops, fontsize=6, zorder=20)
    ax.annotate('', xy=(x_sc, y_sc-1), xytext=(x_sc, y_sc+ly_sc),
                        arrowprops=arrowprops, fontsize=6, zorder=20)
    ax.text(x_sc+lx_sc, y_sc, '$x$', fontsize=9,
                    verticalalignment='bottom', horizontalalignment='center', bbox=my_bbox)
    ax.text(x_sc-0.5, y_sc+ly_sc, '$z$', verticalalignment='center',
                    horizontalalignment='right', fontsize=9, bbox=my_bbox)


for i in range(N):
    for j in range(M):
        axs[i, j].plot([-200, 200], [0, 0], 'k', linewidth=0.5)
        axs[i, j].plot([0, 0], [-200, 200], 'k', linewidth=0.5)
        axs[i, j].text(-0.5, 0, 'O', fontsize=10, horizontalalignment='right',
                    verticalalignment='bottom', bbox=dict(boxstyle="round", fc="w", alpha=0.2))
        draw_sc = 0
        if i == 0 and j == 0:
            x_sc, y_sc = -15, 0
            lx_sc, ly_sc = 5, 5
            draw_ax_sc(axs[i, j],x_sc, y_sc, lx_sc, ly_sc)


x_sc, y_sc = 10, 5
lx_sc, ly_sc = 10, 5
# draw_ax_sc(ax2, x_sc, y_sc, lx_sc, ly_sc, fmt='.0f', shift_y_txt=.5)





#Draw numerical res box
from matplotlib.patches import Circle, Wedge, Polygon, Arc, Rectangle

dS5=5*np.array(sR.dConfig['box'])/np.array(sR.dConfig['res'])*100
pos_rect=[20,10]
rect = Rectangle(pos_rect, -dS5[0], -dS5[2], ec="none",color='k', linewidth=0.5,zorder=1)
# ax2.add_patch(rect)
# ax2.text(pos_rect[0]-dS5[0]/2,10.5,r'5 $\delta x$',horizontalalignment='center')
# ax2.text(pos_rect[0]+.5,8.5,r'5 $\delta z$')

# draw gravity vector
x_sc, y_sc = 25, 12
lx_sc, ly_sc = 10, 5
vecgy = y_sc
vecgx = x_sc
lg = 8
dxT = 0.2*lg
dyT = 0.1 * lg


for j in range(4):
    slopeinrad = runExps[j].dictE['Slope'] * 3.15 / 180
    axs[j, 2].quiver(vecgx, vecgy, lg*np.sin(slopeinrad), -lg*np.cos(slopeinrad),
                    width=0.004, linewidth=1, angles='xy', scale_units='xy', scale=1)
    axs[j, 2].plot([vecgx, vecgx], [vecgy, vecgy-1.1*lg], 'k', linewidth=1)

    arc = Arc([vecgx, vecgy], lg*1.4, lg*1.4, theta1=270,
            theta2=270+np.rad2deg(slopeinrad), linewidth=0.5, color='k')
    axs[j, 2].add_patch(arc)
    axs[j, 2].text(vecgx+lg*slopeinrad+dxT, vecgy-lg *
                np.cos(slopeinrad), r'$\vec{g}$', fontsize=10)
    axs[j, 2].text(vecgx+lg*slopeinrad+dxT, vecgy-lg*np.cos(slopeinrad)+5*dyT,
                r'$i$='+format(runExps[j].dictE['Slope'], '1.0f')+r'$^\circ $', fontsize=10)
plt.pause(0.1)


ax_draw.set_xlim([0, 1])
ax_draw.set_ylim([0, 1])

NsR = len(selectedRuns)
indContrst = 0
ls = ['-.', '-', '-.', '--']
c = [0.9, 1.1, 1., 1.2]
if Nrun > 7:
    shiftExp = -1
else:
    shiftExp=5

Vini = np.zeros((len(selectedRuns)))
# init Vini

for sR, i in zip(selectedRuns, range(len(selectedRuns))):
    nF = int(sR.dConfig['nFrames'])
    ifile=3
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    sR.nYplot=int(sR.resY//2)
    contrs = sR.findContourPhi(level=0.5)
    V = area(contrs[0])
    Vini[i] = V

    k = 0
    for ifile in [6, 12, nF]:
        ax = axs[i, k]
        h1, l1 = [], []
    
        contrs=np.array(dictArt[fignames[i]]['experiment']['free_surface'][str(ifile)])
        
        for cont in contrs:
            contn=np.array(cont)
            [line2D] = ax.plot(smooth(contn[:, 0],10), smooth(contn[:, 1],10), linestyle='-', color='purple', linewidth=1.2, alpha=0.7, label="Exp.")
            

        if ifile < nF:
            # runExps[i].plotField(ax, np.max(
                # [(ifile - 3)* 10 - shiftExp, 0]), vmin=0, vmax=1, cmap=cmap)
            norm=(runExps[i].Ux[(ifile - 3)* 10 - shiftExp]**2+runExps[i].Uy[(ifile - 3)* 10 - shiftExp]**2)**0.5
        # else: 
        #     runExps[i].plotField(ax, -1, vmin=0, vmax=1, cmap=cmap)
        #     norm=(runExps[i].Ux[-1]**2+runExps[i].Uy[-1]**2)**0.5    
            contrs_vel = d6py.Tools.findContours(runExps[i].X, runExps[i].Y, norm.T, 0.01)

            for cont in contrs_vel:
                [line2D_vel]=ax.plot(smooth(cont[:, 0],10)*100, smooth(cont[:, 1],10)*100, linestyle='--', color='purple', linewidth=0.8, alpha=0.8, label="Exp.")


        sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))

        sR.plotContour(ax, levels=[0.5], linewidths=1.2, linestyles='-.')
        
        sR.calculateNormVelocity()
        sR.normV[np.where(sR.normV==0)]=np.nan
        contours=d6py.Tools.findContours(sR.grid_x[:,0, 0], sR.grid_z[0,0,:], sR.normV[:,sR.nYplot,:], 0.01)
        if ifile<nF:
            for cont in contours:
                [line2D_vel_mod] = axs[i, k].plot(smooth(cont[:, 0],5)*100, smooth(cont[:, 1],5)*100, linestyle=':', color='k', linewidth=0.9, alpha=1, label="limit-mod")
            
        V = area(sR.findContourPhi(level=0.5)[0])

        lost = (Vini[i]-V)/Vini[i]*100

        mod = 'velocity'
        im = sR.opyfPointCloudColoredScatter(
            ax, nvec=6000, mute=True, vmin=0, vmax=1, s=0.8, cmap=cmap, rasterized=True, mod=mod)
        sR.plotDoor(ax, alpha=0.5)

        h = sR.CS.legend_elements(str(sR.dimSim)+"D~-~ \mu= " + toS(sR.dConfig['mu'], 2))[0]
        l = [str(sR.dimSim)+r"D Sim. free surf."]
        h1, l1 = h1+h, l1+l

        sR.plotDoor(ax, alpha=0.5)

        k += 1
        
    plt.pause(0.1)        


# axs[1,-1].legend(h1+ [line2D_vel_mod]  + [line2D] + [line2D_vel] , l1+ [r"3D Sim. static-flowing trans."] +
#            [r"Exp. free surface"] + [r"Exp. static-flowing trans."] , fontsize=7, framealpha=0.5, loc=1)


[x, y, X, Y] = axs[0, 0].get_position().bounds
plt.figtext(x+X/2, y+Y+0.025, r'$t=0.2$ s',
            fontsize=9, horizontalalignment='center')
[x, y, X, Y] = axs[0, 1].get_position().bounds
plt.figtext(x+X/2, y+Y+0.025, r'$t=0.6$ s',
            fontsize=9, horizontalalignment='center')
[x, y, X, Y] = axs[0, 2].get_position().bounds
plt.figtext(x+X/2, y+Y+0.025, r'$t_f$',
            fontsize=9, horizontalalignment='center')

fig.legend(h1+ [line2D_vel_mod]  + [line2D] + [line2D_vel] , l1+ [r"3D Sim. static-flowing"] + [r"Exp. free surf."] + [r"Exp. static-flowing"], fontsize=9,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=4,bbox_to_anchor=(0.03, -0.005, 0.4, 0.2))
plt.pause(0.1)
fig.savefig(driveFolder+"/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular-all/dry-granular/doc/article/images/used_images/G.pdf", dpi=150)
# sys.stdout = sys.__stdout__
print(r'\includegraphics{test_savefig_pdf.pdf}')


# %%

