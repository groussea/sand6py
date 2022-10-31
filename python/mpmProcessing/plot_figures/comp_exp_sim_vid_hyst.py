#%%
%matplotlib qt
from template_runout import *
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Arc
import importlib
importlib.reload(d6py)
os.chdir(os.path.dirname(os.path.abspath(__file__)))

Nrun=7
driveFolder = '/media/gauthier/Data-Gauthier/Gauthier'
driveFolder = '/media/gauthier/DataSSD'
paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(
    maind6OutFolder, 4)

R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0, runNumber=Nrun, dimSim=3,  muRigid=0.23, mu=0.44,keyWord='infNorm',fps=150)
scale = 0.01  # 1cm
for sR in R1:
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
OR2 = OR2-np.mean(HE[-int(len(HE)/2):-30])/runExp.dictE['scale']

runExp.video.scaleData(
    metersPerPx=runExp.dictE['scale']/scale,
    framesPerSecond=runExp.dictE['fps'], origin=OR2)
# Interpolate values on field
# xfMss, VxfMss, tfs, maxVss, times, dt= loadRunOutMod(sR,hth=hth,vth=vth)

#%%

plt.close('all')
sc=0.4
fig, axs = plt.subplots(1,1,figsize=(16*sc,9.*sc))

L = (runExp.dictE['L']-0.01) / runExp.scaleLength
H = (runExp.dictE['H']+0.02) / runExp.scaleLength
N = 2
w_fig = 15

Y_lim = [-0.015/runExp.scaleLength, H]
X_lim = [-L-0.2, runExp.xmax * H]
w_axs = 0.85
h_axs = 0.3
w_s = 0.12
w_s2 = 0.02

    
x_sc_ax, y_sc_ax = 30, 5
lx_sc_ax, ly_sc_ax = 5, 5
    

x_sc, y_sc = 47, 12
lx_sc, ly_sc = 10, 5
vecgy = y_sc
vecgx = x_sc
lg = 8
dxT = 0.2*lg
dyT = 0.1 * lg
slopeinrad = runExp.dictE['Slope'] * 3.15 / 180

ideb=int(0.2*sR.fps)

sR.loadVTK(ideb)
sR.calculatePointsLayer() 
ptsI=sR.selectedInd.copy()
ind = np.random.choice(np.arange(len(sR.pointsLayer)), 10000, replace=False)

im=sR.opyfPointCloudColoredScatter(axs,mod='velocity',vmin=0,vmax=1,cmap=cmap, rasterized=True, s=2,nvec=10000,ind=ind)

def calibAxes(axs,L,w_s,N,H,h_axs,w_axs,X_lim,Y_lim):
    axs.plot([-L, -L], [-10, 2 * H], '-k', linewidth=3)
    [x, y, X, Y] = axs.get_position().bounds
    axs.set_position(
        [w_s, 0.4, w_axs, h_axs])
    axs.set_xlim([X_lim[0], 60])
    axs.set_ylim(Y_lim)
    axs.set_aspect('equal', adjustable='box')
    axs.set_yticklabels([])
    axs.set_xticklabels([])
    axs.grid()
    axs.plot([-200, 200], [0, 0], 'k', linewidth=0.5)
    axs.plot([0, 0], [-200, 200], 'k', linewidth=0.5)
    axs.text(-0.5, 0, 'O', fontsize=8, horizontalalignment='right',
                    verticalalignment='bottom', bbox=dict(boxstyle="round", fc="w", alpha=0.2))
    
calibAxes(axs,L,w_s,N,H,h_axs,w_axs,X_lim,Y_lim)

[x, y, X, Y] = axs.get_position().bounds

folderout="/media/gauthier/DataSSD/programs/gitLab/dry-granular/doc/article/video/exp_vs_sim/"+fignames[Nrun]+'_hyst/'
opyf.mkdir2(folderout)
ls2= [':',':',':']
ls = ['-.', '-.', '-.', '--']
c = [1.2, 1.0, 0.8, 1.2]
c2 = [0.9, 1. , 0.6]
col = [ 'k','yellowgreen',  'cornflowerblue']
#%
for i_t in range(ideb,sR.nF):
# for i_t in range(250,251):
#%  
    axs.cla()
    axs.cla()
    h1, l1 = [], [] 
    calibAxes(axs,L,w_s,N,H,h_axs,w_axs,X_lim,Y_lim)
    draw_ax_sc(axs, x_sc_ax, y_sc_ax, lx_sc_ax, ly_sc_ax)    
    
    axs.quiver(vecgx, vecgy, lg*np.sin(slopeinrad), -lg*np.cos(slopeinrad),
                     width=0.004, linewidth=1, angles='xy', scale_units='xy', scale=1)
    axs.plot([vecgx, vecgx], [vecgy, vecgy-1.1*lg], 'k', linewidth=1)

    arc = Arc([vecgx, vecgy], lg*1.4, lg*1.4, theta1=270,
              theta2=270+np.rad2deg(slopeinrad), linewidth=0.5, color='k')
    axs.add_patch(arc)
    axs.text(vecgx+lg*slopeinrad+dxT, vecgy-lg *
                   np.cos(slopeinrad), r'$\vec{g}$', fontsize=7)
    axs.text(vecgx+lg*slopeinrad+dxT, vecgy-lg*np.cos(slopeinrad)+5*dyT,
                   r'$i$='+format(runExp.dictE['Slope'], '1.0f')+r'$^\circ $', fontsize=7)
    


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
        [line2D] = axs.plot(cont[:, 0],  cont[:, 1], linestyle='-', color='purple', linewidth=1.4, alpha=0.7, label="exp.")

    # axs.imshow(vis, extent=runExp.video.paramPlot['extentFrame'])
    # runExp.plotField(axs, np.max(
            # [i_t-30 - shiftExp, 0]), vmin=0, vmax=1, cmap=cmap)

    norm = (runExp.Ux[i_t-30 - shiftExp]**2 +
            runExp.Uy[i_t-30 - shiftExp]**2)**0.5
    if i_t>40:
        contrs_vel_exp = d6py.Tools.findContours(
            runExp.X, runExp.Y, norm.T, 0.01)
        for cont in contrs_vel_exp:
            [line2D_vel] = axs.plot(cont[:, 0]*100, cont[:, 1] *100, linestyle='--', color='purple', linewidth=1, alpha=0.7, label="Exp.")
    s_t=(i_t-ideb)/sR.fps
    
    
    # for i_t in range(ideb,sR.nF):

    txt1=fig.text(0.02,0.52,r'$t='+format(s_t,'1.1f')+'$ s',fontsize=11,)

    for sR, i in zip(R1, range(len(R1))):
        sR.loadVTK(i_t)

        sR.selectedInd=ptsI
        sR.pointsLayer = sR.pointsp[ptsI, :]
        # im=sR.opyfPointCloudColoredScatter(axs,mod='velocity',vmin=0,vmax=1,cmap=cmap, rasterized=True, s=2,nvec=10000,ind=ind)
        if i == 0:
            sR.plotContour(axs, colors=col[i % 4],levels=[0.5], linewidths=c[i % 4], linestyles=ls[i])
        if i == 1:
            sR.plotContour(axs, levels=[0.5],colors=col[i % 4],  linewidths=c[i % 4], linestyles=ls[i],)
        sR.calculateNormVelocity()
        sR.normV[np.where(sR.normV == 0)] = np.nan
        contours = d6py.Tools.findContours(sR.grid_x[:, 0, 0], sR.grid_z[0, 0, :], sR.normV[:, sR.nYplot, :], 0.01)

        for cont in contours:
            linesC= axs.plot(cont[:, 0]*100, cont[:, 1]*100, linestyle=ls2[i], color=col[i % 4], linewidth=c[i % 4], alpha=1, label="limit-mod")

        h = sR.CS.legend_elements()[0]
        l = [r"sim. free surf. "]
        h1, l1 = h1+h, l1+l
        l = [r"sim. stat.-flow. $\Delta_{\mu_{hyst}}=" + format(sR.dConfig['delta_mu_start'],'0.1f') + r"$ | $I_*=5 \times 10^{-3}$"]
        h1, l1 = h1+ linesC, l1+l 

    l_new=[l1[0],l1[1],l1[2],l1[3],r"exp. free surf.", r"exp. stat.-flow."]

    h_new=[h1[0],h1[1],h1[2],h1[3],line2D,line2D_vel]

    dy=0.05
    ly=0.25
    legend1=fig.legend(h_new[0:2]  , l_new[0:2] , fontsize=10,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=2,bbox_to_anchor=(0.15,ly, 0.5, 0.3))
    legend2=fig.legend(h_new[2:4]  , l_new[2:4] , fontsize=10,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=2,bbox_to_anchor=(0.15, ly-dy, 0.5, 0.3))
    legend3=fig.legend(h_new[4:6]  , l_new[4:6] , fontsize=10,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=2,bbox_to_anchor=(0.15, ly-2*dy, 0.5, 0.3))

    
    sR.plotDoor(axs)
    sR.plotDoor(axs)
    # axs.legend(h_new_2, l_new_2, fontsize=7, framealpha=1, loc=9)
    # axs.legend(h_new_1, l_new_1, fontsize=7, framealpha=1, loc=9)

    # fig.show()
    [x, y, X, Y] = axs.get_position().bounds
    txt2=fig.text(x+X/2, y+Y+0.1, 'B15 - hysteresis - 3D simulations versus experiment',
            fontsize=12, horizontalalignment='center')
    fig.savefig(folderout+format(i_t,'03.0f')+".png",dpi=300)
    txt1.remove()
    txt2.remove()
    legend1.remove()
    legend2.remove()
    legend3.remove()
    
    # plt.pause(0.5)



