
#%%
%matplotlib qt
from template_runout import *

driveFolder = '/media/gauthier/Data-Gauthier/Gauthier'
driveFolder = '/media/gauthier/DataSSD'
# maind6OutFolder = '/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/'
paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(
    maind6OutFolder, 4)

#%
Nrun=7
run=str(Nrun)

fignames = ['G00', 'G05', 'G10', 'G15', 'B00', 'B05', 'B10', 'B15', 'B20']
run = fignames[Nrun]
# R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0., runNumber=Nrun, dimSim=3, delta_mu_start=0., mu=0.44, keyWord='try2')

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth  
plt.close('all')

# fig, ax, ax12= generate_fig()
# fig2, ax2, ax22= generate_fig()




mainExpFolder = driveFolder + \
        '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'
cmap=plt.get_cmap('inferno')
ic=0
npl=1
hth,vth=0.007, 0.001 # Threshold on runout height in (m) and on final time en m/s
tfsM,tfsE=[],[]
outf='/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/transient_runouts_B_runs/'

# for Nrun in range(7,7+npl):
Nrun=7
run=fignames[Nrun]
# R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0, runNumber=Nrun, dimSim=3, delta_mu_start=0.0, muRigid=0.23, mu=0.44,keyWord='W_8')
R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0, runNumber=Nrun, dimSim=3, delta_mu_start=0.0, muRigid=0.23, mu=0.44,keyWord='infNorm')
sR=R1[1]
# Interpolate values on field
# xfMss, VxfMss, tfs, maxVss, times, dt= loadRunOutMod(sR,hth=hth,vth=vth)
sR.nF = int(sR.dConfig['nFrames'])
sR.fps = int(sR.dConfig['fps'])
meanIs, stdIs,Iqlows, Iq50s, Iqhighs, histIs =[],[], [], [], [], []
jj=0
ideb=int(0.2*sR.fps)
for i in range(ideb,sR.nF):
# for i in range(8,9):
    sR.loadVTK(i)
    # sR.calc_I()
    # sR.calculateNormVelocity()
    # sR.IFieldsup0=[np.where(sR.IField>0.02)]
    # normV2D=np.flipud(sR.normV[:,sR.nYplot,:].T)
    # I_moving=sR.IField[np.where((normV2D>0.01)*(sR.IField>1e-5))]
    norm=np.linalg.norm(sR.vel[:, :],axis=1)
    # I_moving=sR.inertia[np.where((norm>0.001))]
    # I_moving=sR.inertia[np.where((norm>0.01)*(sR.inertia>5e-3))]
    I_moving=sR.inertia[np.where((sR.inertia>5e-4))]
    # I_moving=sR.IField[np.where((normV2D>0.01)*(sR.IMgamma>1))]
    # I_moving=sR.IField[np.where((sR.IField>1e-6))]

    # I_moving=sR.IField[np.where((sR.IField>1e-4)*(sR.IMgamma>10)*(sR.P[:,sR.nYplot,:]>1e-3))]
    if len(I_moving)<10:
        I_moving=0
    meanI=np.mean(I_moving)
    meanIs.append(meanI)
    stdI=np.std(I_moving)
    stdIs.append(stdI)
    Iqlow=np.quantile(I_moving,0.3)
    Iqlows.append(Iqlow)
    Iq50=np.quantile(I_moving,0.5)
    Iq50s.append(Iq50)
    Iqhigh=np.quantile(I_moving,0.7)
    Iqhighs.append(Iqhigh)
    histI=np.histogram(I_moving, bins=30, range=[0,0.1], normed=None, weights=None, density=True)
    # ax.plot(histI[0])
    histIs.append(histI[0])
    t=(i-3)/sR.fps
    # if jj==0:
    #     ax.plot(t,Iqlow,'ro',label=' I quantile 40')
    #     ax.plot(t,Iq50,'bd',label=' I median')
    #     ax.plot(t,Iqhigh,'g+',label=' I quantile 60')
    #     ax.plot(t,meanI,'k*',label=' I averaged')
    # else:
    #     ax.plot(t,Iqlow,'ro')
    #     ax.plot(t,Iq50,'bd')
    #     ax.plot(t,Iqhigh,'g+')
    #     ax.plot(t,meanI,'k*')
    jj+=1
    # ax.plot([t,t],[meanI-stdI, meanI+stdI],'k+')
    # ax.plot([t,t],[Iq20, Iq80],'b+')
    # ax.plot(t,meanI-stdI,'k+')

#%%
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Arc
def my_formatter(x, pos):
    """Format 1 as 1, 0 as 0, and all values whose absolute values is between
    0 and 1 without the leading "0." (e.g., 0.7 is formatted as .7 and -0.4 is
    formatted as -.4)."""
    val_str = '{:g}'.format(x)
    if np.abs(x) > 0 and np.abs(x) < 1:
        return val_str.replace("0", "", 1)
    else:
        return val_str
major_formatter = FuncFormatter(my_formatter)

plt.close('all')
fig, axs = plt.subplots(1,2,figsize=(7,2.))
times=np.array([(i-ideb)/sR.fps for i in range(ideb,sR.nF)])
histIs=np.array(histIs)
# axs[1].imshow(np.flipud(histIs.T),vmin=0, vmax=100, extent=[times[0],times[-1],0,0.1],cmap='hot_r',aspect='auto',interpolation='gaussian')
# ax.plot(times,Iqlows,'r--',label=' I quantile 40', lw=1)

axs[1].plot(times,Iq50s,'k--',label=' $I$ median', lw=1)
axs[1].fill_between(times,Iqlows, Iqhighs, color='blue', alpha=0.5,label=' $I$  [30,\,70]-quantiles ')
# ax.plot(times,Iqhighs,'g--',label=' I quantile 60', lw=1)
# ax.plot(times,meanIs,'k-',label=' I averaged')
s_t=1.1
i_t=np.where(times>s_t)[0][0]

sR.loadVTK(i_t+ideb)
sR.scLength(0.11)
sR.calc_I()
im=sR.opyfPointCloudColoredScatter(axs[0],mod='inertia',vmin=0,vmax=0.3,alpha=0.9,cmap='magma_r', rasterized=True, s=2,nvec=10000)
cbaxes = fig.add_axes([0.5, 0.9, 0.48, 0.05])

cb = fig.colorbar(im, cax=cbaxes, orientation='horizontal', extend='both')
cb.set_label('$I$', fontsize=10)

#draw experimental profile on top

mainExpFolder = driveFolder + \
    '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'
runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
scale=0.11
runExp1.scLength(scale)
# runExp1.plotDepthProfile(axs[0],i_t,c='purple',lw=1)    


# [line2D] = axs[0].plot(smooth(contn[:, 0],1), smooth(contn[:, 1],1), linestyle='-', color='purple', linewidth=1.5, alpha=0.7, label="Exp.")
#%

# axs[0].axis('equal')
ypos=0.18
H0=11
axs[0].set_position([0.5, ypos, 0.48, 0.55])
axs[0].set_xlim([-22/H0,90/H0])
axs[0].set_ylim([-1/H0,9/H0])
axs[1].set_position([0.1, ypos, 0.30, 0.75])
axs[1].set_xlabel('$t$ [s]')
axs[1].legend()
axs[1].set_ylabel('$I$')
axs[1].set_xlim([0,times[-1]])
axs[1].set_ylim([0.001,0.086])
axs[0].set_xlabel(r'$x/H_0$')
axs[0].set_ylabel(r'$y/H_0$')
axs[0].grid()
axs[1].grid()
plt.show()
axs[0].text(24/H0,6./H0,r'$t_{I_{70,max}}='+format(s_t,'1.1f')+'$ s',fontsize=14,)
axs[1].plot([1.11],[0.052], 'ko',lw=1,ms=4)
axs[1].plot([1.11,1.4],[0.052,0.052], 'k-',lw=1)
axs[1].text(1.45,0.052,r'$I_{70,max}$',fontsize=12,verticalalignment='center')

x_sc, y_sc = 7, 0.7
lx_sc, ly_sc = 2, 1
vecgy = y_sc
vecgx = x_sc
lg = 0.5
dxT = 0.2*lg
dyT = 0.1 * lg
slopeinrad = 15 * 3.15 / 180
axs[0].quiver(vecgx, vecgy, lg*np.sin(slopeinrad), -lg*np.cos(slopeinrad),
                width=0.005, linewidth=1, angles='xy', scale_units='xy', scale=1)
axs[0].plot([vecgx, vecgx], [vecgy, vecgy-1.1*lg], 'k', linewidth=1)

arc = Arc([vecgx, vecgy], lg*1.4, lg*1.4, theta1=270,
        theta2=270+np.rad2deg(slopeinrad), linewidth=0.5, color='k')
axs[0].add_patch(arc)
axs[0].text(vecgx+lg*slopeinrad+dxT, vecgy-lg *
            np.cos(slopeinrad), r'$\vec{g}$', fontsize=10)
axs[0].text(vecgx+lg*slopeinrad+dxT, vecgy-lg*np.cos(slopeinrad)+5*dyT,
            r'$i$='+format(runExp1.dictE['Slope'], '1.0f')+r'$^\circ $', fontsize=10)


for a in axs:
    a.xaxis.set_major_formatter(major_formatter)
    a.yaxis.set_major_formatter(major_formatter)

fig.text(0.04,0.9,'(a)')
fig.text(0.45,0.9,'(b)')
fig.savefig("/media/gauthier/DataSSD/programs/gitLab/dry-granular/doc/article/figures/evolution_I.pdf",dpi=300)

#%%





#%%

# plt.subplots()


# plt.imshow(sR.normV[:,sR.nYplot,:].T)
# plt.figure()
# plt.imshow(sR.IField)

# sR.plotI(ax)

    # tfsM.append(tfs[0])
    # VxfMssP=np.insert(VxfMss[0][:-1:],0,0)
    # timesP=np.insert(times[:-2:]+dt/2,0,0)
    # plot=ax.plot(timesP, VxfMssP,'-',c=cmap(ic/npl+0.05),lw=0.9, label=run)


    # opyf.mkdir2(outf+run)