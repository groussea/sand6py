
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
Nrun = 8
scale = 0.01  # 1cm
runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
runExp1.scLength(scale)
mu = runExp1.dictE['mu']



delta_mu=0.

# R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=delta_mu, fps=15, nSamples=2, I0_start=0.00, mu=0.38, delta_mu_start=0., runNumber=Nrun, dimSim=2,substeps=80,viscosity=0.0)
R1 = d6py.NumericalRun('/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/Run_08_3D_Door_mu=0.48_muRigid=0.18_W_3.0cm_glass-beads-0.47mm_Slope=20deg_delta_mu=0.26_substeps_80_fracH=0.8_I0_start=0.0000_delta_mu_start=0.00test-scaling/')
R1=[R1]
R2 = d6py.NumericalRun('/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/Run_08_3D_Door_mu=0.48_muRigid=0.18_W_6.0cm_glass-beads-0.47mm_Slope=20deg_delta_mu=0.26_substeps_80_fracH=0.8_I0_start=0.0000_delta_mu_start=0.00test-scaling/')
R2=[R2]

# R2, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=delta_mu, mu=0.48, fps=15, nSamples=2, I0_start=0.00, delta_mu_start=0., runNumber=Nrun, dimSim=2, substeps=80, viscosity=0)





R3, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=delta_mu, muRigid=0.18, mu=0.48, fps=15, nSamples=2, runNumber=Nrun, dimSim=3)

R4, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=delta_mu, fps=15, nSamples=2, I0_start=0.00, mu=0.48, delta_mu_start=0., runNumber=Nrun, dimSim=2,substeps=80,viscosity=0.0)

# delta_mu=0.22

# R4, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=delta_mu, muRigid=0.0, mu=0.48, fps=15, nSamples=2, runNumber=Nrun, dimSim=2)

# delta_mu=0.22
# R3, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=delta_mu, muRigid=0.18, mu=0.38, fps=15, nSamples=2,  runNumber=Nrun, dimSim=3)

# R4, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=delta_mu, muRigid=0.18, mu=0.48, fps=15, nSamples=2, runNumber=Nrun, dimSim=3)

selectedRuns=[R1[0],R2[0],R3[0],R4[0]]
# selectedRuns = [Ref[0], selectedRuns[0]]


for sR in selectedRuns:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])
    sR.nYplot=3

#%%

runExp1.loadVideo('/media/gauthier/Samsung_T5/MPM_data/Collapse_experiment/Video_src', mute=True,display=False)

# # plt.close('all')
step, shift, Ntot = 8, 5, (int(runExp1.dictE['nFrames'])+20)*2
runExp1.video.set_vecTime(starting_frame=int(runExp1.dictE['framedeb']) - 100, step=step, shift=shift, Ntot=Ntot)

OR2 = runExp1.dictE['OR2']


BW, vecXE, HE = mask_collapses(opyf.Tools.convertToGrayScale(runExp1.video.cropFrameInit), OR2, 0.8, runExp1.dictE['scale'])
OR2[1] = OR2[1]-np.mean(HE[-int(len(HE)/2):-30])/runExp1.dictE['scale']


runExp1.video.scaleData( metersPerPx=runExp1.dictE['scale'] / scale, framesPerSecond=runExp1.dictE['fps'], origin=OR2)


#%%


from matplotlib.ticker import FuncFormatter
major_formatter = FuncFormatter(d6py.Tools.my_formatter)

plt.ion()
plt.rcParams['font.family'] = 'serif'

plt.rc('text', usetex=True)

plt.close('all')

cmapg = opyf.custom_cmap.make_cmap_customized(Palette='green')

colors = [(33./255, 66./255, 99./255, 0.05),
          (1, 1, 0.3, 0.8), (0.8, 0, 0, 0.8), (0, 0, 0, 0.8)]
position = [0, 0.1, 0.5, 1]
cmap = opyf.make_cmap(colors, position)

# cmap.set_over(color='g')
cmap.set_under(alpha=0)
L = 0.22/0.01
H =0.15/0.01

w_fig = 13.5
N = 3
M = 4
fig, axs = plt.subplots(N, M, dpi=142, figsize=(w_fig * 0.39, 15 * 0.39))

# fig = plt.figure(dpi=142, figsize=(w_fig * 0.39, 11 * 0.39))

Y_lim = [0, H]
X_lim = [-L, L]
w_axs = 0.195
h_axs = 0.135
w_s = 0.07
w_s2 = 0.04
for i in range(N):
    for j in range(M):
        axs[i, j].plot([-L-1, -L-1], [-2, 2 * H], '-k', linewidth=4)
        [x, y, X, Y] = axs[i, j].get_position().bounds
        axs[i, j].set_position(
            [w_s + (w_axs + w_s2)*j, 0.18+(h_axs + 0.055)*(N - i), w_axs, h_axs])
        axs[i, j].xaxis.set_major_formatter(major_formatter)
        axs[i, j].yaxis.set_major_formatter(major_formatter)
        if i == 0:
            axs[i, j].set_xlim([-23, 20])
            axs[i, j].set_ylim([-1.5, 12])
        if i == 1:
            axs[i, j].set_xlim([-23, 22])
            axs[i, j].set_ylim([-1.5, 12])
        if i == 2:
            axs[i, j].set_xlim([-23, 60])
            axs[i, j].set_ylim([-1.5, 12])
        # axs[i, j].set_aspect('equal', adjustable='box')
        if j>0:
            axs[i, j].set_yticklabels([])
        if j == 0:
            axs[i, j].set_ylabel('z (cm)',fontsize=7)
        
        if i == 2:
            axs[i, j].set_xlabel('x (cm)', fontsize=7)
            axs[i, j].xaxis.set_label_coords(0.5,-0.15)
        # axs[i, j].set_xticklabels([])

# draw lines to separate experimental plots
ax_draw = fig.add_axes([0, 0, 1., 1.], zorder=-10)
ax_draw.set_axis_off()
ax_draw.grid()
[x, y, X, Y] = axs[2, 1].get_position().bounds
X_line = x+X+w_s2/2

[x, y, X, Y] = axs[2, 3].get_position().bounds
Y_sep = y-0.03
ax_draw.plot([X_line, X_line], [Y_sep, 0.98], linewidth=0.5, color='k')
# ax_draw.plot([0.02, 0.98], [Y_sep, Y_sep], linewidth=0.5, color='k')
# Add a big axe to draw the final state entirely
ax2 = fig.add_axes([w_s, 0.07, x + X - w_s, Y_sep-0.2], zorder=-10)

ax2.set_xlim([-L, 60])
ax2.plot([-L, -L], [-10, 2 * H], '-k', linewidth=2)
ax2.set_ylim([-1.5,11])
ax2.set_ylabel('z (cm)', fontsize=8)
ax2.set_xlabel('x (cm)', fontsize=8)
ax2.xaxis.set_major_formatter(major_formatter)
ax2.yaxis.set_major_formatter(major_formatter)
# ax2.set_aspect('equal', adjustable='box')
# draw scale and axes
ax2.yaxis.set_label_coords(-0.04,0.4)

for i in range(N):
    for j in range(M):
        axs[i, j].plot([-200, 200], [0, 0], 'k', linewidth=0.5)
        axs[i, j].plot([0, 0], [-200, 200], 'k', linewidth=0.5)
        axs[i, j].text(-0.5, 0, 'O', fontsize=8, horizontalalignment='right',
                       verticalalignment='bottom', bbox=dict(boxstyle="round", fc="w", alpha=0.2))


my_bbox = dict(fc="w", alpha=0.3)
ax2.plot([-200, 200], [0, 0], 'k', linewidth=0.5)
ax2.plot([0, 0], [-200, 200], 'k', linewidth=0.5)
ax2.text(-0.5, 0, 'O', fontsize=7, horizontalalignment='right',
         verticalalignment='bottom', bbox=my_bbox)

x_sc, y_sc = 60, 2
lx_sc, ly_sc = 10, 2


# #Draw numerical res box
from matplotlib.patches import Circle, Wedge, Polygon, Arc, Rectangle


x_sc, y_sc = 4, 11
lx_sc, ly_sc = 10, 5
vecgy = y_sc
vecgx = x_sc
lg = 6
dxT = 0.2*lg
dyT = 0.1 * lg
slopeinrad = runExp1.dictE['Slope'] * 3.15 / 180

for j in [0]:
    d6py.Tools.pltGravity(axs[0, j],x_sc,y_sc,lx_sc,ly_sc,vecgx,vecgy,lg,dxT,dyT,slopeinrad)


# write Experiments and Model
[x, y, X, Y] = axs[0, 0].get_position().bounds
plt.figtext(x+X+w_s2/2, y+Y+0.07, '2D',
            fontsize=10, horizontalalignment='center')

[x, y, X, Y] = axs[0, 2].get_position().bounds
plt.figtext(x+X+w_s2/2, y+Y+0.07, r'3D - $\mu_w=0.18$',
            fontsize=10, horizontalalignment='center')

# for j in range(0,4):
#     [x, y, X, Y] = axs[0, j].get_position().bounds
#     plt.figtext(x + X / 2, y + Y + 0.03, '$\mu_1$='+ toS(selectedRuns[j].dConfig['mu'], 2) , fontsize=8, horizontalalignment='center')
    # plt.figtext(x + X / 2+0.06, y + Y + 0.03, '$\mu_w$='+ toS(selectedRuns[j].dConfig['muRigid'], 2) , fontsize=8, horizontalalignment='center')



ax_draw.set_xlim([0, 1])
ax_draw.set_ylim([0, 1])


plt.show()

# %
k = 0
NsR = len(selectedRuns)
indContrst = 0
ls = ['--', '-.', '--', '-.']
c = [0.9, 1., 0.8, 1.2]
SR = selectedRuns[0:4]
if Nrun > 7:
    shiftExp = -1
else:
    shiftExp = 5
    
Vini = np.zeros((len(SR)))



# init Vini
ifile=3
for sR, i in zip(SR, range(len(SR))):
    ax = axs[k, i]
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    contrs = sR.findContourPhi(level=0.5)
    V = area(contrs[0])
    Vini[i] = V

for ifile in [ 6, 9 , 15]:
    
    ax = axs[k, 0]
    h1, l1 = [], []
    time = (ifile-3) / sR.dConfig['fps'] 
    # write times

    [x, y, X, Y] = axs[k, 0].get_position().bounds
    plt.figtext(0.02, y + Y *1.05, 't='+format(time,'1.1f') +' s', fontsize=8, fontweight='bold')


    for sR, i in zip(SR, range(len(SR))):

        ax = axs[k, i]
        sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
        [line2D]=runExp1.plotDepthProfile(ax, np.max(
            [(ifile-3) * 10 - 5, 0]), linestyle='-', color='purple', linewidth=0.9, label="Experience",zorder=1)
        sR.plotContour(ax, levels=[0.5], linewidths=c[i % 4], linestyles=ls[i % 4], colors=[ cmapg((NsR - i) / (NsR + indContrst))],zorder=2)
        
        V = area(sR.findContourPhi(level=0.5)[0])

        lost=(Vini[i]-V)/Vini[i]*100

        [x, y, X, Y] = axs[k, i].get_position().bounds
        plt.figtext(x+X*0.55, y + 1.05*Y, r'$\epsilon$='+format(lost,'1.1f') +r'\%',color=(0.33,0,0), fontsize=7)
        
        im = sR.opyfPointCloudColoredScatter(
            ax, nvec=8000, mute=True, vmin=0, vmax=1, s=0.03, cmap=cmap, rasterized=True)

        sR.plotDoor(ax, alpha=0.5)
        if k ==0 and (i==0 or i==2):
            h1, l1 = [[line2D], ["Exp."]]
            h = sR.CS.legend_elements()[0]   
            l = [r"$\mu_1$= " + toS(sR.dConfig['mu'],2)]
            h1, l1 = h1+h, l1+l
        if k == 0 and (i == 1 or i == 3):
            h = sR.CS.legend_elements()[0]   
            l = [r"$\mu_1$= " + toS(sR.dConfig['mu'],2)]
            h1, l1 = h1+h, l1+l           
            fig.legend(h1 , l1, fontsize=7,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=3,bbox_to_anchor=(x-X-1.5*w_s2, 0.9, 0.4, 0.2))




    ax = axs[k, 0]
    # sR.plotDoor(ax, alpha=0.5)
    # plt.pause(0.2)
    k += 1


ifile = nF


    

runExp1.video.readFrame(np.max([int(runExp1.dictE['framedeb']) - 100 + (ifile-3) * (
    100) - shiftExp*10, int(runExp1.dictE['framedeb']) - 100]))
vis = opyf.Render.CLAHEbrightness(runExp1.video.vis, 150)
ax2.imshow(vis, extent=runExp1.video.paramPlot['extentFrame'])

sys.path.append('/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python/imageProcessing/')

BW=mask_collapses2(runExp1.video.vis,0.5)


mat = np.array(BW, dtype=np.float32)
runExp1.video.set_gridToInterpolateOn(stepGrid=1)
x = runExp1.video.vecX
y=np.flipud(runExp1.video.vecY)
contrs = d6py.Tools.findContours(x, y, np.flipud(mat).T, 10)



[line2D] = ax2.plot(contrs[0][:,0], contrs[0][:,1],linestyle='-', color='purple', linewidth=0.9, alpha=0.5,label="Exp.")

# runExp1.plotDepthProfile(ax2, np.max(
#             [(ifile - 3)* 10 - 5, 0]), linestyle='-', color='k', linewidth=0.8, label="Exp.")

[fantom]=ax.plot(1,1,'-',alpha=0)     
h1, l1 = [[line2D], ["Exp."]]
h1, l1 = h1+[fantom], l1+['']
for sR, i in zip(selectedRuns, range(len(selectedRuns))):
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    sR.plotContour(ax2, levels=[0.5], linewidths=c[i % 4], linestyles=ls[i % 4], colors=[
                   cmapg((NsR-i)/(NsR+indContrst))])
    h = sR.CS.legend_elements()[0]
    
    l = [str(sR.dimSim) + r"D~-~$\mu_1$= " + toS(sR.dConfig['mu'],2)]
    h1, l1 = h1+h, l1+l
ax2.legend(h1 , l1, fontsize=7, framealpha=0.5, loc=1,ncol=3)
time = (ifile-3) / sR.dConfig['fps'] 
[x, y, X, Y] = ax2.get_position().bounds
plt.figtext(0.02, y + 1.05*Y , 't='+format(time,'1.1f') +' s - Final state', fontsize=8,fontweight='bold')
[x, y, X, Y] = axs[2, 2].get_position().bounds
cbaxes = fig.add_axes([0.18, y-0.09, 0.70, 0.02])

cb = fig.colorbar(im, cax=cbaxes, orientation='horizontal', extend='both')
cb.set_label('Velocity norm [m/s]', fontsize=8)


plt.show()
plt.pause(0.2)
# fig.savefig("test_savefig_pdf_5_mm.pdf", dpi=300)
fig.savefig("/media/gauthier/Data-Gauthier/Gauthier/TAF/TAF_inria/INRIA_current_work/GitLab/dry-granular-all/dry-granular/doc/article/images/used_images/Run_08_3D-2D.pdf", dpi=150)
# sys.stdout = sys.__stdout__
print(r'\includegraphics{test_savefig_pdf.pdf}')



# %%



