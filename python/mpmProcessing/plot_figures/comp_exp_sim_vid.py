#%%
%matplotlib qt
from template_runout import *
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Arc
import importlib
importlib.reload(d6py)
os.chdir(os.path.dirname(os.path.abspath(__file__)))

Nrun=5
driveFolder = '/media/gauthier/Data-Gauthier/Gauthier'
driveFolder = '/media/gauthier/DataSSD'
paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(
    maind6OutFolder, 4)

R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0, runNumber=Nrun, dimSim=3, delta_mu_start=0.0, muRigid=0.23, mu=0.44,keyWord='infNorm',fps=150)
sR=R1[0]

scale = 0.01  # 1cm
sR.scLength(scale)
mainExpFolder = driveFolder + \
    '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'
runExp = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)


    

vidFolder='/media/gauthier/DataSSD/temp/Video_src/'
runExp.loadVideo(vidFolder, mute=True)
step, shift, Ntot = 8, 5, (int(runExp.dictE['nFrames'])+50)*2
runExp.video.set_vecTime(starting_frame=int(
    runExp.dictE['framedeb']) - 100, step=step, shift=shift, Ntot=Ntot)
OR2 = runExp.dictE['OR2']
runExp.scLength(scale)


BW, vecXE, HE = mask_collapses(opyf.Tools.convertToGrayScale(
    runExp.video.cropFrameInit), OR2, 0.8, runExp.dictE['scale'])
OR2[1] = OR2[1]-np.mean(HE[-int(len(HE)/2):-30])/runExp.dictE['scale']

runExp.video.scaleData(
    metersPerPx=runExp.dictE['scale']/scale,
    framesPerSecond=runExp.dictE['fps'], origin=OR2)
# Interpolate values on field
# xfMss, VxfMss, tfs, maxVss, times, dt= loadRunOutMod(sR,hth=hth,vth=vth)

#%%

plt.close('all')
sc=0.4
fig, axs = plt.subplots(2,1,figsize=(16*sc,9.*sc))

L = (runExp.dictE['L']-0.01) / runExp.scaleLength
H = (runExp.dictE['H']+0.02) / runExp.scaleLength
N = 2
w_fig = 15

Y_lim = [-0.015/runExp.scaleLength, H]
X_lim = [-L-0.2, runExp.xmax * H]
w_axs = 0.75
h_axs = 0.3
w_s = (1-w_axs)/2
w_s2 = 0.02

    
x_sc_ax, y_sc_ax = 23, 5
lx_sc_ax, ly_sc_ax = 5, 5
    

x_sc, y_sc = 35, 12
lx_sc, ly_sc = 10, 5
vecgy = y_sc
vecgx = x_sc
lg = 8
dxT = 0.15*lg
dyT = 0.1 * lg
slopeinrad = runExp.dictE['Slope'] * 3.15 / 180

ideb=int(0.2*sR.fps)

sR.loadVTK(ideb)
sR.calculatePointsLayer() 
ptsI=sR.selectedInd.copy()
ind = np.random.choice(np.arange(len(sR.pointsLayer)), 10000, replace=False)
im=sR.opyfPointCloudColoredScatter(axs[0],mod='velocity',vmin=0,vmax=1,cmap=cmap, rasterized=True, s=2,nvec=10000,ind=ind)
cbaxes = fig.add_axes([0.1, 0.1, 0.8, 0.03])

cb = fig.colorbar(im, cax=cbaxes, orientation='horizontal', extend='both')
cb.set_label('norm velocity [m/s]', fontsize=10)

def calibAxes(axs,L,w_s,N,H,h_axs,w_axs,X_lim,Y_lim):
    for i in range(N):
        axs[i].plot([-L, -L], [-10, 2 * H], '-k', linewidth=3)
        [x, y, X, Y] = axs[i].get_position().bounds
        axs[i].set_position(
            [w_s, 0.17+(h_axs +0.13) * (N-i-1), w_axs, h_axs])
        axs[i].set_xlim([X_lim[0], 45])
        axs[i].set_ylim(Y_lim)
        axs[i].set_aspect('equal', adjustable='box')
        axs[i].set_yticklabels([])
        axs[i].set_xticklabels([])
        axs[i].grid()
        axs[i].plot([-200, 200], [0, 0], 'k', linewidth=0.5)
        axs[i].plot([0, 0], [-200, 200], 'k', linewidth=0.5)
        axs[i].text(-0.5, 0, 'O', fontsize=8, horizontalalignment='right',
                        verticalalignment='bottom', bbox=dict(boxstyle="round", fc="w", alpha=0.2))
        
calibAxes(axs,L,w_s,N,H,h_axs,w_axs,X_lim,Y_lim)

[x, y, X, Y] = axs[1].get_position().bounds
plt.figtext(x+X/2, y+Y+0.045, 'Experiment',
            fontsize=9, horizontalalignment='center')
[x, y, X, Y] = axs[0].get_position().bounds
plt.figtext(x+X/2, y+Y+0.045, r'3D simulation - $\mu= '+format(sR.dConfig['mu'])+'$',
            fontsize=9, horizontalalignment='center')

# fig.text(0.02, y+Y+0.045,'G15 - 2.7 mm granules - 15°',fontsize=11)
fig.text(0.02, y+Y+0.045,'B05 - 0.5 mm beads - 5°',fontsize=11)
folderout="/media/gauthier/DataSSD/programs/gitLab/dry-granular/doc/article/video/exp_vs_sim/"+fignames[Nrun]+'/'
opyf.mkdir2(folderout)
for i_t in range(ideb,sR.nF):
# for i_t in range(30,31):
#%
    axs[0].cla()
    axs[1].cla()
    calibAxes(axs,L,w_s,N,H,h_axs,w_axs,X_lim,Y_lim)
    draw_ax_sc(axs[0], x_sc_ax, y_sc_ax, lx_sc_ax, ly_sc_ax)    
    
    axs[0].quiver(vecgx, vecgy, lg*np.sin(slopeinrad), -lg*np.cos(slopeinrad),
                     width=0.004, linewidth=1, angles='xy', scale_units='xy', scale=1)
    axs[0].plot([vecgx, vecgx], [vecgy, vecgy-1.1*lg], 'k', linewidth=1)

    arc = Arc([vecgx, vecgy], lg*1.4, lg*1.4, theta1=270,
              theta2=270+np.rad2deg(slopeinrad), linewidth=0.5, color='k')
    axs[0].add_patch(arc)
    axs[0].text(vecgx+lg*slopeinrad+dxT, vecgy-lg *
                   np.cos(slopeinrad), r'$\vec{g}$', fontsize=7)
    axs[0].text(vecgx+lg*slopeinrad+dxT, vecgy-lg*np.cos(slopeinrad)+5*dyT,
                   r'$i$='+format(runExp.dictE['Slope'], '1.0f')+r'$^\circ $', fontsize=7)
    
    
    sR.loadVTK(i_t)
    sR.selectedInd=ptsI
    sR.pointsLayer = sR.pointsp[ptsI, :]
    im=sR.opyfPointCloudColoredScatter(axs[0],mod='velocity',vmin=0,vmax=1,cmap=cmap, rasterized=True, s=2,nvec=10000,ind=ind)
    sR.plotContour(axs[0], levels=[0.5], linewidths=1.4, linestyles='-.')
    sR.calculateNormVelocity()
    sR.normV[np.where(sR.normV == 0)] = np.nan
    contours = d6py.Tools.findContours(sR.grid_x[:, 0, 0], sR.grid_z[0, 0, :], sR.normV[:, sR.nYplot, :], 0.01)

    for cont in contours:
        [line2D_vel_mod] = axs[0].plot(cont[:, 0]*100, cont[:, 1] * 100,  linestyle=':', color='k', linewidth=1.1, alpha=1, label="limit-mod")



    shiftExp=4
    runExp.video.readFrame(np.max([int(runExp.dictE['framedeb']) - 100 + i_t*10-300 - shiftExp*10, int(runExp.dictE['framedeb']) - 100]))
    vis = opyf.Render.CLAHEbrightness(runExp.video.vis, 100)
    BW = mask_collapses2(runExp.video.vis, 0.8)

    mat = np.array(BW, dtype=np.float32)
    runExp.video.set_gridToInterpolateOn(stepGrid=1)
    x = runExp.video.vecX
    y = np.flipud(runExp.video.vecY)
    contrs = d6py.Tools.findContours(x, y, np.flipud(mat).T, 10)
    for cont in contrs:
        [line2D] = axs[1].plot(cont[:, 0],  cont[:, 1], linestyle='-', color='purple', linewidth=1.4, alpha=0.7, label="exp.")

    axs[1].imshow(vis, extent=runExp.video.paramPlot['extentFrame'])
    runExp.plotField(axs[1], np.max(
            [i_t-30 - shiftExp, 0]), vmin=0, vmax=1, cmap=cmap)

    norm = (runExp.Ux[i_t-30 - shiftExp]**2 +
            runExp.Uy[i_t-30 - shiftExp]**2)**0.5
    if i_t>40:
        contrs_vel_exp = d6py.Tools.findContours(
            runExp.X, runExp.Y, norm.T, 0.01)
        for cont in contrs_vel_exp:
            [line2D_vel] = axs[1].plot(cont[:, 0]*100, cont[:, 1] *100, linestyle='--', color='purple', linewidth=1, alpha=0.7, label="Exp.")
    s_t=(i_t-ideb)/sR.fps
    sR.plotDoor(axs[0])
    sR.plotDoor(axs[1])
    txt1=fig.text(0.02,0.73,r'$t='+format(s_t,'1.1f')+'$ s',fontsize=11,)
    txt2=fig.text(0.02,0.3,r'$t='+format(s_t,'1.1f')+'$ s',fontsize=11,)
    # for i_t in range(ideb,sR.nF):
    l_new_2 = [r"sim. free surf.", r"sim. stat.-flow."]
    h = sR.CS.legend_elements(str(sR.dimSim)+"D~-~\mu_{RB}=" + toS(
    sR.dConfig['muRigid'], 2) + " ~-~ \mu= " + toS(sR.dConfig['mu'], 2) + " - \phi")[0]
    h_new_2 = [h[0], line2D_vel_mod]

    l_new_1 = [r"exp. free surf.", r"exp. stat.-flow."]
    h_new_1 = [ line2D,line2D_vel]

    axs[0].legend(h_new_2, l_new_2, fontsize=7, framealpha=1, loc=9)
    axs[1].legend(h_new_1, l_new_1, fontsize=7, framealpha=1, loc=9)

    # fig.show()

    fig.savefig(folderout+format(i_t,'03.0f')+".png",dpi=300)
    
    txt1.remove()
    txt2.remove()
    # plt.pause(0.5)



