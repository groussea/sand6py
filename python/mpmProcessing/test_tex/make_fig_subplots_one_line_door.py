
# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Author : Gauthier Rousseau
import opyf  # from opyflow library some rendering function may be employed
import d6py
import sys
sys.path.append(
    '/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python/imageProcessing')

from Tools_collapses import mask_collapses, mask_collapses2
import matplotlib.pyplot as plt
from d6py.Tools import *
import argparse
# plt.ioff()
# intialize exteral packages
import sys
import os
import numpy as np
import matplotlib

%matplotlib qt5
# plt.ioff()
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth

# sys.path.append('/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python')
# sys.path.append('/media/gauthier/
# Data-Gauthier/programs/gitHub/opyflow')

outDictFolder=    "/media/gauthier/Data-Gauthier/Gauthier/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular-all/dry-granular/data/outputs_experiments_and_mpm/"
fignames=['G00', 'G05', 'G10', 'G15', 'B00', 'B05', 'B10', 'B15', 'B20']


# init_article_dictionary(mainExpFolder,outDictFolder)     
JSONpath = outDictFolder + "article_dict.json"
in_file = open(JSONpath,"r")
dictArt= json.load(in_file) 
in_file.close()

# sys.stdout = open(os.devnull, 'w')
# intialize font type and size
plt.rcParams['font.size'] = 8.0
plt.rcParams['xtick.labelsize'] = 8.0
plt.rcParams['ytick.labelsize'] = 8.0
plt.rcParams['ytick.labelsize'] = 8.0
plt.rcParams['axes.linewidth'] = 1
# print(r'\includegraphics{test_2.pdf}')
driveFolder = '/media/gauthier/Data-Gauthier/Gauthier'
maind6OutFolder = '/media/gauthier/Samsung_T51/sand6_sorties/sand6_out/'
# maind6OutFolder = '/home/gauthier/sorties_sand6/'
paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(
    maind6OutFolder, 4)
# %matplotlib qt5
parser = argparse.ArgumentParser()
parser.add_argument('--ny', help="number of elements in y direction",default=6)
parser.add_argument('--path', help="out path", default="")
# args = parser.parse_args()

Nrun=7
scale = 0.01  # 1cm
mainExpFolder = driveFolder + \
    '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'
runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
runExp1.scLength(scale)
mu = runExp1.dictE['mu']


muw=0.
# 0.0],[0.1,0.15]
i0=0.005
if Nrun < 4:
    mu=0.75
else:
    mu=0.44
#%%

R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, mu=mu, delta_mu=0., runNumber=Nrun, dimSim=3, delta_mu_start=0, keyWord='try2')

R2, selectedDict = d6py.whereSand6OutFromParms(listNumRun, mu=mu, delta_mu=0., runNumber=Nrun, dimSim=3, delta_mu_start=0, keyWord='door')

#%%
# selectedRuns=[Ref[-1],selectedRuns[-1],selectedRuns2[-1]]
selectedRuns = [R1[0], R2[0]] #8
# selectedRuns = [R1[0], R2[0]] #7

for sR in selectedRuns:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])


runExp1.loadVideo(
    '/media/gauthier/Data-Gauthier/Gauthier/TAF/TAF_inria/MPM-data/Collapse_Experiment/Video_src/', mute=True)
plt.close('all')
step, shift, Ntot = 8, 5, (int(runExp1.dictE['nFrames'])+20)*2
runExp1.video.set_vecTime(starting_frame=int(
    runExp1.dictE['framedeb']) - 100, step=step, shift=shift, Ntot=Ntot)
OR2 = runExp1.dictE['OR2']


BW, vecXE, HE = mask_collapses(opyf.Tools.convertToGrayScale(
    runExp1.video.cropFrameInit), OR2, 0.8, runExp1.dictE['scale'])
OR2[1] = OR2[1]-np.mean(HE[-int(len(HE)/2):-30])/runExp1.dictE['scale']


runExp1.video.scaleData(
    metersPerPx=runExp1.dictE['scale']/scale, framesPerSecond=runExp1.dictE['fps'], origin=OR2)
dictArt[fignames[Nrun]]['experiment']['OR2b']=OR2
f = open(JSONpath, "w")
json.dump(dictArt, f, indent=4)
f.close() 
# %%

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

# cmap.set_over(color='g')
cmap.set_under(alpha=0)
L = runExp1.dictE['L'] / runExp1.scaleLength-1
H = (runExp1.dictE['H']+0.02 ) / runExp1.scaleLength

w_fig = (L + runExp1.xmax *
        runExp1.dictE['H'] / runExp1.scaleLength) / 83.4 * 8

w_fig = 16
N = 2  # lignes
M = 3 # colonnes
fig, axs = plt.subplots(N, M, dpi=142, figsize=(w_fig * 0.5, 1.3))
Y_lim = [-0.015/runExp1.scaleLength, H]
X_lim = [-L-0.2, runExp1.xmax * H]
w_axs = 0.23
h_axs = 0.6
w_s = 0.02
w_s2 = 0.01




for i in range(N):
    for j in range(M):
        axs[i, j].plot([-L, -L], [-10, 2 * H], '-k', linewidth=2)
        axs[i, j].set_anchor('SW')
        [x, y, X, Y] = axs[i, j].get_position().bounds
        axs[i, j].set_position(
            [w_s+(w_axs + w_s2) * j, (h_axs+0.05 ) * (N-i-1)+0.25, w_axs*3, h_axs])
        axs[i, j].set_xlim([X_lim[0], 7.5])
        axs[i, j].set_ylim(Y_lim)
        axs[i, j].set_aspect('equal')
        axs[i, j].set_yticklabels([])
        axs[i, j].set_xticklabels([])



# draw lines to separate experimental plots
ax_draw = fig.add_axes([0, 0, 1, 1], zorder=-10)
ax_draw.grid()
ax_draw.set_axis_off()
X_line = w_s + (w_axs + 0.01)
for i in range(N): 
    axs[i,-1].set_xlim([X_lim[0], 50])
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
    ax.text(x_sc-1, y_sc+ly_sc-1, '$z$', verticalalignment='center',
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
slopeinrad = runExp1.dictE['Slope'] * 3.15 / 180

for j in [2]:
    axs[1, j].quiver(vecgx, vecgy, lg*np.sin(slopeinrad), -lg*np.cos(slopeinrad),
                    width=0.004, linewidth=1, angles='xy', scale_units='xy', scale=1)
    axs[1, j].plot([vecgx, vecgx], [vecgy, vecgy-1.1*lg], 'k', linewidth=1)

    arc = Arc([vecgx, vecgy], lg*1.4, lg*1.4, theta1=270,
            theta2=270+np.rad2deg(slopeinrad), linewidth=0.5, color='k')
    axs[1, j].add_patch(arc)
    axs[1, j].text(vecgx+lg*slopeinrad+dxT, vecgy-lg *
                np.cos(slopeinrad), r'$\vec{g}$', fontsize=10)
    axs[1, j].text(vecgx+lg*slopeinrad+dxT, vecgy-lg*np.cos(slopeinrad)+5*dyT,
                r'$i$='+format(runExp1.dictE['Slope'], '1.0f')+r'$^\circ $', fontsize=10)

#%
# write times

# plt.figtext(x+X/2, y+Y+0.045, 't=0 s',
            # fontsize=9, horizontalalignment='center')


# for j in range(1,3):
#     [x, y, X, Y] = axs[0, j].get_position().bounds
#     plt.figtext(x + X / 2-0.06, y + Y +0.015, '$\mu_1$='+ toS(selectedRuns[j-1].dConfig['mu'], 2) , fontsize=8, horizontalalignment='center')
#     if j==2:
#         plt.figtext(x + X / 2+0.06, y + Y +0.015, '$\mu_w$='+ toS(selectedRuns[j-1].dConfig['muRigid'], 2) , fontsize=8, horizontalalignment='center')

ax_draw.set_xlim([0, 1])
ax_draw.set_ylim([0, 1])


arrowprops = dict(
    arrowstyle="-",
    color='black')


axs[1, 0].annotate('uplifting gate', xy=(2, 15), xytext=(-3, 14.5),
                arrowprops=arrowprops, fontsize=6, zorder=20, verticalalignment='center', horizontalalignment='right', bbox=my_bbox)

arrowprops = dict(
    arrowstyle="->",
    color='black')

axs[1, 0].annotate('', xy=(4, 15), xytext=(4, 10),
                arrowprops=arrowprops, fontsize=6, zorder=20, verticalalignment='center', horizontalalignment='right')



# %
k = 0
NsR = len(selectedRuns)
indContrst = 0
ls = ['-', '-.', '-.', '--']
ls2= ['--',':']
c = [0.9, 1.1, 1., 1.2]
col = ['black', 'blue']
SR = selectedRuns[0:2]
if Nrun > 7:
    shiftExp = -1
else:
    shiftExp=5

Vini = np.zeros((len(SR)))
# init Vini
ifile=3
for sR, i in zip(SR, range(len(SR))):
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    sR.nYplot=int(sR.resY//2)
    contrs = sR.findContourPhi(level=0.5)
    V = area(contrs[0])
    Vini[i] = V

if 'free_surface' not in dictArt[fignames[Nrun]]['experiment'].keys():
    dictArt[fignames[Nrun]]['experiment']['free_surface']={}

for ifile in [6, 12, nF]:
    ax = axs[0, k]
    h1, l1 = [], []
    if ifile<nF:   
        runExp1.video.readFrame(np.max([int(runExp1.dictE['framedeb']) - 100 + (ifile-3) * (100) - shiftExp*10, int(runExp1.dictE['framedeb']) - 100]))
    else:
        runExp1.video.readFrame(runExp1.video.vec[-1])   
        
    vis = opyf.Render.CLAHEbrightness(runExp1.video.vis, 150)
    BW=mask_collapses2(runExp1.video.vis,0.6)

    mat = np.array(BW, dtype=np.float32)
    runExp1.video.set_gridToInterpolateOn(stepGrid=1)
    x = runExp1.video.vecX
    y=np.flipud(runExp1.video.vecY)
    contrs = d6py.Tools.findContours(x, y, np.flipud(mat).T, 10)
    # save contour in the json file
    
    dictArt[fignames[Nrun]]['experiment']['free_surface'][str(ifile)]=[]

    # contrs=np.array(dictArt[fignames[Nrun]]['experiment']['free_surface'][str(ifile)])
    
    for cont in contrs:
        dictArt[fignames[Nrun]]['experiment']['free_surface'][str(ifile)].append(cont.tolist())
        contn=np.array(cont)
        # [line2D] = ax.plot(smooth(contn[:, 0],10), smooth(contn[:, 1],10), linestyle='-', color='purple', linewidth=1.5, alpha=0.7, label="Exp.")
        

    # ax.imshow(vis, extent=runExp1.video.paramPlot['extentFrame'])
    # if ifile < nF:
        # runExp1.plotField(ax, np.max(
        #     [(ifile - 3)* 10 - shiftExp, 0]), vmin=0, vmax=1, cmap=cmap)
        # norm=(runExp1.Ux[(ifile - 3)* 10 - shiftExp]**2+runExp1.Uy[(ifile - 3)* 10 - shiftExp]**2)**0.5
    # else: 
    #     runExp1.plotField(ax, -1, vmin=0, vmax=1, cmap=cmap)
    #     norm=(runExp1.Ux[-1]**2+runExp1.Uy[-1]**2)**0.5    
        # contrs_vel = d6py.Tools.findContours(runExp1.X, runExp1.Y, norm.T, 0.01)

        # for cont in contrs_vel:
            # [line2D_vel]=ax.plot(smooth(cont[:, 0],10)*100, smooth(cont[:, 1],10)*100, linestyle='--', color='purple', linewidth=0.8, alpha=0.7, label="Exp.")
            # axs[1, k].plot(smooth(cont[:, 0],10)*100, smooth(cont[:, 1],10)*100, linestyle='--', color='purple', linewidth=0.8, alpha=0.7, label="Exp.")

    for sR, i in zip(SR, range(len(SR))):
        axt = axs[1, k]
        sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))

        sR.plotContour(axt, levels=[0.5], linewidths=c[i % 4], linestyles=ls[i % 4], colors=col[i % 4], alpha=0.6+i*0.4)
        sR.calculateNormVelocity()
        sR.normV[np.where(sR.normV==0)]=np.nan
        contours=d6py.Tools.findContours(sR.grid_x[:,0, 0], sR.grid_z[0,0,:], sR.normV[:,sR.nYplot,:], 0.01)
        if ifile<nF:
            for ii in [0,1]:
                for cont in contours:
                    [line2D_vel_mod] = axs[ii, k].plot(smooth(cont[:, 0],5)*100, smooth(cont[:, 1],5)*100, linestyle=ls2[i % 4], color=col[i % 4], linewidth=c[i % 4], alpha=0.6+i*0.4, label="limit-mod")
                    
            
        V = area(sR.findContourPhi(level=0.5)[0])

        lost = (Vini[i]-V)/Vini[i]*100

        mod = 'velocity'
        # if i>0:
            # im = sR.opyfPointCloudColoredScatter( axt, nvec=6000, mute=True, vmin=0, vmax=1, s=0.8, cmap=cmap, rasterized=True, mod=mod)
        sR.plotDoor(axt, alpha=0.5)
        # for cont in contrs:
        #     contn=np.array(cont)
        #     [line2D] = axt.plot(contn[:, 0], contn[:, 1], linestyle='-', color='purple', linewidth=1.5, alpha=0.7, label="Exp.")

        axt.set_yticklabels([])
        axt.set_xticklabels([])

        h = sR.CS.legend_elements(str(sR.dimSim)+"D~-~ \mu= " + toS(sR.dConfig['mu'], 2))[0]
        if i<1:
            l = [r"Sim. free surf. - $\mu_D=0.18$"]
        else:
            l = [r"Sim. free surf. - $\mu_D=1.0$"]
        h1, l1 = h1+h, l1+l
        if i<1:
            l = [r"Sim. static-flowing trans. - $\mu_D=0.18$"]
        else:
            l = [r"Sim. static-flowing trans. - $\mu_D=0.18$"]
        h1, l1 = h1+ [line2D_vel_mod], l1+l

    sR.plotDoor(ax, alpha=0.5)
    k += 1
    
axs[0,0].remove()
axs[0,1].remove()
axs[0,2].remove()

# axs[1,-1].legend(h1+ [line2D_vel_mod]  + [line2D] + [line2D_vel] , l1+ [r"3D Sim. static-flowing trans."] +
#            [r"Exp. free surface"] + [r"Exp. static-flowing trans."] , fontsize=7, framealpha=0.5, loc=1)

fig.legend(h1  , l1 , fontsize=6,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=4,bbox_to_anchor=(0.03, 0.005, 0.4, 0.2))

[x, y, X, Y] = axs[1, 0].get_position().bounds
# plt.figtext(x+X/2, y+Y+0.045, r'$t=0.2$ s',
#             fontsize=9, horizontalalignment='center')
# [x, y, X, Y] = axs[1, 1].get_position().bounds
# plt.figtext(x+X/2, y+Y+0.045, r'$t=0.6$ s',
#             fontsize=9, horizontalalignment='center')
# [x, y, X, Y] = axs[1, 2].get_position().bounds
# plt.figtext(x+X/2, y+Y+0.045, r'$t_f$',
#             fontsize=9, horizontalalignment='center')

f = open(JSONpath, "w")
json.dump(dictArt, f, indent=4)
f.close() 
fig.set_size_inches()
fig.savefig("/media/gauthier/Data-Gauthier/Gauthier/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular-all/dry-granular/doc/article/images/used_images/"+fignames[Nrun]+"_door.pdf", dpi=150)
# sys.stdout = sys.__stdout__
print(r'\includegraphics{test_savefig_pdf.pdf}')


# %%

