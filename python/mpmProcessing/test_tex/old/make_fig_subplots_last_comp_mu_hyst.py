
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

# intialize exteral packages
import sys
import os
import numpy as np
import matplotlib
# matplotlib.use("Qt5Agg")
# sys.path.append('/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python')
# sys.path.append('/media/gauthier/Data-Gauthier/programs/gitHub/opyflow')

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
paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(
    maind6OutFolder, 4)
%matplotlib qt5
Nrun = 6

scale = 0.01  # 1cm

mainExpFolder = driveFolder + \
    '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'
runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
runExp1.scLength(scale)
mu = runExp1.dictE['mu']

if Nrun <4:
    dmu, dmus = 0.15, 0.1
else:
    dmu, dmus = 0.2, 0.15
muw=0.
# 0.0],[0.1,0.15]
i0=0.005
mu=0.65
dmu=-0.05


R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, muRigid=0.18,mu=0.44, delta_mu=0., runNumber=Nrun, dimSim=3, delta_mu_start=0)

R2, selectedDict = d6py.whereSand6OutFromParms(listNumRun, muRigid=0.18, mu=0.54, delta_mu=0., runNumber=Nrun, dimSim=3, delta_mu_start=0.1)

#%%
# selectedRuns=[Ref[-1],selectedRuns[-1],selectedRuns2[-1]]
# selectedRuns = [R1[0], R2[0]] #8
selectedRuns = [R2[0], R2[-1]] #7

for sR in selectedRuns:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])

runExp1.loadVideo(
    '/media/gauthier/Samsung_T51/MPM_data/Collapse_experiment/Video_src', mute=True)
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
# %%
plt.close('all')
plt.ion()
plt.rcParams['font.family'] = 'serif'
plt.rc('text', usetex=True)
plt.close('all')

cmapg = opyf.custom_cmap.make_cmap_customized(Palette='green')

colors = [(33./255, 66./255, 99./255, 0.05),
          (1, 1, 0.3, 0.9), (0.8, 0, 0, 0.9), (0, 0, 0, 0.9)]
position = [0, 0.1, 0.5, 1]
cmap = opyf.make_cmap(colors, position)

# cmap.set_over(color='g')
cmap.set_under(alpha=0)
L = runExp1.dictE['L'] / runExp1.scaleLength
H = (runExp1.dictE['H']+0.02 ) / runExp1.scaleLength

w_fig = (L + runExp1.xmax *
         runExp1.dictE['H'] / runExp1.scaleLength) / 83.4 * 8

w_fig = 15
N = 3
M = 3
fig, axs = plt.subplots(N, M, dpi=142, figsize=(w_fig * 0.39, 13 * 0.39))

# fig = plt.figure(dpi=142, figsize=(w_fig * 0.39, 11 * 0.39))

Y_lim = [-0.015/runExp1.scaleLength, H]
X_lim = [-L-0.2, runExp1.xmax * H]
w_axs = 0.28
h_axs = 0.2
w_s = 0.09
w_s2 = 0.02
for i in range(N):
    for j in range(M):

        axs[i, j].plot([-L, -L], [-10, 2 * H], '-k', linewidth=2)
        [x, y, X, Y] = axs[i, j].get_position().bounds
        axs[i, j].set_position(
            [w_s+(w_axs + w_s2) * j, 0.2+(h_axs - 0.03) * (N-i), w_axs, h_axs])
        axs[i, j].set_xlim([X_lim[0], 7.5])
        axs[i, j].set_ylim(Y_lim)
        axs[i, j].set_aspect('equal', adjustable='box')
        axs[i, j].set_yticklabels([])
        axs[i, j].set_xticklabels([])

# draw lines to separate experimental plots
ax_draw = fig.add_axes([0, 0, 1, 1], zorder=-10)
ax_draw.grid()
ax_draw.set_axis_off()
X_line = w_s + (w_axs + 0.01)
[x, y, X, Y] = axs[2, 2].get_position().bounds
Y_sep = y-0.03
ax_draw.plot([X_line, X_line], [Y_sep, 0.98], linewidth=0.5, color='k')
ax_draw.plot([0.02, 0.98], [Y_sep, Y_sep], linewidth=0.5, color='k')

# Add a big axe to draw the final state entirely
ax2 = fig.add_axes([w_s, 0.00, x + X - w_s, Y_sep-0.03], zorder=-10)
ax2.set_yticklabels([])
ax2.set_xticklabels([])
ax2.set_xlim([X_lim[0], 49])
ax2.plot([-L, -L], [-10, 2 * H], '-k', linewidth=2)
ax2.set_ylim(Y_lim)
ax2.set_aspect('equal', adjustable='box')
# draw scale and axes

my_bbox = dict(fc="w", alpha=0.3)

def draw_ax_sc(ax,x_sc, y_sc,lx_sc, ly_sc,fmt='.0f',shift_y_txt=0.5):
    ax.plot([x_sc, x_sc+lx_sc],
                    [y_sc, y_sc], '-+k', linewidth=1)
    ax.plot([x_sc, x_sc], [y_sc, y_sc + ly_sc],
                    '-+k', linewidth=1)
    ax.text(x_sc + lx_sc / 2, y_sc + shift_y_txt, format(lx_sc, fmt) + ' cm', fontsize=6,
                    horizontalalignment='center', bbox=my_bbox)
    ax.text(x_sc - 1, y_sc + ly_sc / 2, format(ly_sc, fmt) + ' cm', fontsize=6,
                    horizontalalignment='right', verticalalignment='center', bbox=my_bbox)
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


for i in range(M):
    for j in range(M):
        axs[i, j].plot([-200, 200], [0, 0], 'k', linewidth=0.5)
        axs[i, j].plot([0, 0], [-200, 200], 'k', linewidth=0.5)
        axs[i, j].text(-0.5, 0, 'O', fontsize=8, horizontalalignment='right',
                       verticalalignment='bottom', bbox=dict(boxstyle="round", fc="w", alpha=0.2))
        draw_sc = 0
        if i == 0 and j == 0:
            x_sc, y_sc = -15, 0
            lx_sc, ly_sc = 5, 5
            draw_ax_sc(axs[i, j],x_sc, y_sc,lx_sc, ly_sc)


ax2.plot([-200, 200], [0, 0], 'k', linewidth=0.5)
ax2.plot([0, 0], [-200, 200], 'k', linewidth=0.5)
ax2.text(-0.5, 0, 'O', fontsize=8, horizontalalignment='right',
         verticalalignment='bottom', bbox=my_bbox)

x_sc, y_sc = 0, 5
lx_sc, ly_sc = 10, 5
draw_ax_sc(ax2, x_sc, y_sc, lx_sc, ly_sc, fmt='.0f', shift_y_txt=.5)

#Draw numerical res box
from matplotlib.patches import Circle, Wedge, Polygon, Arc, Rectangle

dS5=5*np.array(sR.dConfig['box'])/np.array(sR.dConfig['res'])*100
pos_rect=[10,10]

rect = Rectangle(pos_rect, -dS5[0], -dS5[2], ec="none",color='k', linewidth=0.5,zorder=1)
ax2.add_patch(rect)
ax2.text(pos_rect[0]-dS5[0]/2,10.5,r'5 $\delta x$',horizontalalignment='center')
ax2.text(pos_rect[0]+0.5,8.5,r'5 $\delta z$')

# draw gravity vector
x_sc, y_sc = -20, 9
lx_sc, ly_sc = 10, 5
vecgy = y_sc
vecgx = x_sc
lg = 8
dxT = 0.2*lg
dyT = 0.1 * lg
slopeinrad = runExp1.dictE['Slope'] * 3.15 / 180
for j in [1]:
    axs[0, j].quiver(vecgx, vecgy, lg*np.sin(slopeinrad), -lg*np.cos(slopeinrad),
                     width=0.008, linewidth=1, angles='xy', scale_units='xy', scale=1)
    axs[0, j].plot([vecgx, vecgx], [vecgy, vecgy-1.1*lg], 'k', linewidth=1)

    arc = Arc([vecgx, vecgy], lg*1.4, lg*1.4, theta1=270,
              theta2=270+np.rad2deg(slopeinrad), linewidth=0.5, color='k')
    axs[0, j].add_patch(arc)
    axs[0, j].text(vecgx+lg*slopeinrad+dxT, vecgy-lg *
                   np.cos(slopeinrad), r'$\vec{g}$', fontsize=7)
    axs[0, j].text(vecgx+lg*slopeinrad+dxT, vecgy-lg*np.cos(slopeinrad)+5*dyT,
                   r'$i$='+format(runExp1.dictE['Slope'], '1.0f')+r'$^\circ $', fontsize=7)


# write times
[x, y, X, Y] = axs[0, 0].get_position().bounds
plt.figtext(0.01, y+Y/2, 't=0 s', fontsize=7)
[x, y, X, Y] = axs[1, 0].get_position().bounds
plt.figtext(0.01, y + Y / 2, 't=0.2 s', fontsize=7)
[x, y, X, Y] = axs[2, 0].get_position().bounds
plt.figtext(0.01, y + Y / 2, 't=0.6 s', fontsize=7)
[x, y, X, Y] = ax2.get_position().bounds
plt.figtext(0.01, y+Y/2, 't=1.6 s', fontsize=7)



#  write Experiments and Model
[x, y, X, Y] = axs[0, 0].get_position().bounds
plt.figtext(x+X/2, y+Y+0.07, 'Experiment',
            fontsize=9, horizontalalignment='center')
[x, y, X, Y] = axs[0, 1].get_position().bounds
plt.figtext(x+X+w_s2/2, y+Y+0.07, '3D simulations',
            fontsize=9, horizontalalignment='center')
plt.figtext(x+X/2, y+Y+0.045, r'Constant $\mu$',
            fontsize=9, horizontalalignment='center')
[x, y, X, Y] = axs[0, 2].get_position().bounds
plt.figtext(x+X/2, y+Y+0.045, r'Hysteretic $\mu$',
            fontsize=9, horizontalalignment='center')

for j in range(1,3):
    [x, y, X, Y] = axs[0, j].get_position().bounds
    plt.figtext(x + X / 2-0.06, y + Y +0.015, '$\mu_1$='+ toS(selectedRuns[j-1].dConfig['mu']-selectedRuns[j-1].dConfig['delta_mu_start'], 2) , fontsize=8, horizontalalignment='center')
    plt.figtext(x + X / 2+0.06, y + Y +0.015, '$\mu_w$='+ toS(selectedRuns[j-1].dConfig['muRigid'], 2) , fontsize=8, horizontalalignment='center')

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


plt.show()

# %
k = 0
NsR = len(selectedRuns)
indContrst = 0
ls = ['--', '-', '-.', '--']
c = [0.9, 1.1, 1., 0.8]
SR = selectedRuns[0:2]
if Nrun > 7:
    shiftExp = -1
else:
    shiftExp=5

# initialisation
ifile=3
    
Vini = np.zeros((len(SR)))

for sR, i in zip(SR, range(len(SR))):

    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    sR.nYplot=int(sR.resY//2)
    contrs = sR.findContourPhi(level=0.5)
    V = area(contrs[0])
    Vini[i] = V

for ifile in [3, 6, 12]:

    ax = axs[k, 0]
    h1, l1 = [], []

    runExp1.video.readFrame(np.max([int(runExp1.dictE['framedeb']) - 100 + (ifile-3) * (
        100) - shiftExp*10, int(runExp1.dictE['framedeb']) - 100]))
    vis = opyf.Render.CLAHEbrightness(runExp1.video.vis, 150)
    BW=mask_collapses2(runExp1.video.vis,0.5)

    mat = np.array(BW, dtype=np.float32)
    runExp1.video.set_gridToInterpolateOn(stepGrid=1)
    x = runExp1.video.vecX
    y=np.flipud(runExp1.video.vecY)
    contrs = d6py.Tools.findContours(x, y, np.flipud(mat).T, 10)
    for cont in contrs:
        [line2D] = ax.plot(cont[:, 0], cont[:, 1], linestyle='-', color='purple', linewidth=1.5, alpha=0.7, label="Exp.")

    ax.imshow(vis, extent=runExp1.video.paramPlot['extentFrame'])
    runExp1.plotField(ax, np.max(
        [(ifile - 3)* 10 - shiftExp, 0]), vmin=0, vmax=1, cmap=cmap)
        

    # runExp1.plotDepthProfile(ax, np.max(
        # [(ifile-3) * 10 - shiftExp, 0]), linestyle='-.', color='k', linewidth=1, label="Experience")


    for sR, i in zip(SR, range(len(SR))):
        ax = axs[k, i + 1]
        sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))

        sR.plotContour(ax, levels=[0.5], linewidths=c[i % 4], linestyles=ls[i % 4], colors=[
                       cmapg((NsR-i)/(NsR+indContrst))])
        mod='velocity'
        im = sR.opyfPointCloudColoredScatter(
            ax, nvec=6000, mute=True, vmin=0, vmax=1, s=0.8, cmap=cmap, rasterized=True, mod=mod)
        # im = sR.plotMu(ax,cmap='viridis',interpolation='gaussian',vmin=0.1, vmax=0.5)
        sR.plotDoor(ax, alpha=0.5)
        for cont in contrs:
            [line2D] = ax.plot(cont[:, 0], cont[:, 1], linestyle='-', color='purple', linewidth=1.5, alpha=0.7, label="Exp.")
        # runExp1.plotDepthProfile(ax, np.max(

        #     [(ifile-3)*10-shiftExp, 0]), linestyle='-.', color='k', linewidth=1, label="Experience")
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        if ifile == nF:
            sR.plotContour(ax2, levels=[0.01], linewidths=c[i % 4], linestyles=ls[i % 4], colors=[
                           cmapg((NsR-i)/(NsR+indContrst))])

        h = sR.CS.legend_elements(str(sR.dimSim)+"D-~\mu_{RB}=" + toS(
            sR.dConfig['muRigid'], 2) + " ~-~ \mu= " + toS(sR.dConfig['mu'], 2) + " - \phi")[0]
        l = [str(sR.dimSim)+r"$D-~\mu_{RB}$=" + toS(
            sR.dConfig['muRigid'], 2) + r"$ ~-~ \mu$= " + toS(sR.dConfig['mu'], 2)]
        h1, l1 = h1+h, l1+l

    ax = axs[k, 0]
    sR.plotDoor(ax, alpha=0.5)
    # plt.pause(0.2)

    k += 1

ifile = nF
h1, l1 = [], []
runExp1.video.vec[-1]
runExp1.video.readFrame(runExp1.video.vec[-1])
vis = opyf.Render.CLAHEbrightness(runExp1.video.vis, 140)
ax2.imshow(vis, extent=runExp1.video.paramPlot['extentFrame'])
# runExp1.plotField(ax2, np.max([ifile * 10 - 5, 0]), vmin=0, vmax=1, cmap=cmap)
# runExp1.plotDepthProfile(ax2, np.max(
#     [(ifile-3) * 10 - 5, 0]), linestyle='-.', color='k', linewidth=1, label="Experience")

# BW=mask_collapses2(runExp1.video.vis,0.5)

if Nrun==8: 
    BW=mask_collapses2(runExp1.video.vis,0.45) #8
else:
    BW=mask_collapses2(runExp1.video.vis,0.5) #8


mat = np.array(BW, dtype=np.float32)
runExp1.video.set_gridToInterpolateOn(stepGrid=1)
x = runExp1.video.vecX
y=np.flipud(runExp1.video.vecY)
contrs = d6py.Tools.findContours(x, y, np.flipud(mat).T, 10)

[line2D] = ax2.plot(contrs[0][:, 0], contrs[0][:, 1], linestyle='-', color='purple', linewidth=1.5, alpha=0.7, label="Exp.")


for sR, i in zip(selectedRuns, range(len(selectedRuns))):
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    V = area(sR.findContourPhi(level=0.5)[0])

    lost=(Vini[i]-V)/Vini[i]*100
    sR.plotContour(ax2, levels=[0.5], linewidths=c[i % 4], linestyles=ls[i % 4], colors=[
                   cmapg((NsR-i)/(NsR+indContrst))])
    h = sR.CS.legend_elements(str(sR.dimSim)+"D-~\mu_{RB}=" + toS(
        sR.dConfig['muRigid'], 2) + " ~-~ \mu= " + toS(sR.dConfig['mu'], 2) + " - \phi")[0]
    if i==0:
        l = [str(sR.dimSim) + r"D ~-~$  \mu=\mu_1$= " + toS(sR.dConfig['mu'], 2)+r'$~\epsilon$='+format(lost,'1.1f') +r'\%']
    if i == 1:
        l = [str(sR.dimSim) + r"D ~-~$\mu_0$= " + toS(sR.dConfig['mu'], 2)+"~-~$\mu_1$= " + toS(sR.dConfig['mu']-sR.dConfig['delta_mu_start'], 2)+r'$~\epsilon$='+format(lost,'1.1f') +r'\%']
    h1, l1 = h1 + h, l1 + l
        
ax2.legend(h1 + [line2D], l1 +[r"Experiment"], fontsize=7, framealpha=0.5, loc=1)


[x, y, X, Y] = axs[2, 2].get_position().bounds
cbaxes = fig.add_axes([0.18, y-0.08, 0.70, 0.03])

cb = fig.colorbar(im, cax=cbaxes, orientation='horizontal', extend='both')
# cb.set_label('Velocity norm [m/s]', fontsize=8)
cb.set_label('Velocity [ m/s ]', fontsize=8)
# cb.set_label('$\mu(I)$ [ - ]', fontsize=8)


[x, y, X, Y] = ax2.get_position().bounds
cbaxes = ax2.set_position([x, y - 0.03, X, Y])
plt.pause(2)
# plt.show()
# fig.savefig("test_savefig_pdf_5_mm.pdf", dpi=300)
fig.savefig("/media/gauthier/Data-Gauthier/Gauthier/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular-all/dry-granular/doc/article/images/used_images/Run_07_hystH.pdf", dpi=150)
# sys.stdout = sys.__stdout__
print(r'\includegraphics{test_savefig_pdf.pdf}')


# %%
