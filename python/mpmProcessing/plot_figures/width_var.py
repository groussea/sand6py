#%%
%matplotlib qt5
from default_2_lines_plot import *


#%%

Nrun=7
scale = 0.01  # 1cm
mainExpFolder = driveFolder + \
    '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'
runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
runExp1.scLength(scale)


if Nrun < 4:
    mu=0.75
else:
    mu=0.44

paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(
    maind6OutFolder, 4)

# R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, mu=mu,  runNumber=Nrun,  keyWord='width')

R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, runNumber=Nrun, keyWord='width_var',mu=0.44)

R2, selectedDict = d6py.whereSand6OutFromParms(listNumRun, runNumber=Nrun, keyWord='W_8',mu=0.44,delta_mu_start=0.00)


#%%
selectedRuns = [R2[0], R1[0], R1[3]] 



for sR in selectedRuns:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])
    print(sR.dConfig['box'][1])

nF=26
    
vidFolder='/media/gauthier/DataSSD/temp/Video_src/'
# vidFolder = '/media/gauthier/Data-Gauthier/Gauthier/TAF/TAF_inria/MPM-data/Collapse_Experiment/Video_src/'
runExp1.loadVideo(vidFolder, mute=True)
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


#%%
plt.close('all')

fig, axs = figure_2lines_template()
draw_gravity_2_lines(axs, runExp1.dictE['Slope'],fontsize_g)

k = 0
NsR = len(selectedRuns)
indContrst = 0
ls = ['-', '-.', '-.', '--']
ls2= [':',':',':']
c = [1.2, 1, 0.8, 1.2]
c2 = [0.9, 1. , 0.6]
col = [ 'k', 'cornflowerblue', 'darkred']
SR = selectedRuns[0:3]
if Nrun > 7:
    shiftExp = -1
else:
    shiftExp=5
lines=[]
Vini = np.zeros((len(SR)))
# init Vini
ifile=3

for sR, i in zip(SR, range(len(SR))):
    sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
    sR.nYplot=int(sR.resY//2)
    contrs = sR.findContourPhi(level=0.5)
    V = area(contrs[0])
    Vini[i] = V


for ifile in [6, 12, nF]:
    print(ifile)
    ax = axs[0, k]
    
    if ifile<nF:   
        h1, l1 = [], [] 
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

    
    # for cont in contrs:
    #     contn=np.array(cont)
    #     [line2D] = axs[1, k].plot(smooth(contn[:, 0],10), smooth(contn[:, 1],10), linestyle='-', color='purple', linewidth=1.5, alpha=0.7, label="Exp.")
        

    if ifile < nF:

        norm=(runExp1.Ux[(ifile - 3)* 10 - shiftExp]**2+runExp1.Uy[(ifile - 3)* 10 - shiftExp]**2)**0.5
  
        contrs_vel = d6py.Tools.findContours(runExp1.X, runExp1.Y, norm.T, 0.01)

        # for cont in contrs_vel:

        #     [line2D_vel]=axs[1, k].plot(smooth(cont[:, 0],10)*100, smooth(cont[:, 1],10)*100, linestyle='--', color='purple', linewidth=0.8, alpha=0.7, label="Exp.")

    for sR, i in zip(SR, range(len(SR))):
        axt = axs[1, k]
        sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
        if i == 0:
            sR.nYplot=3
            sR.plotContour(axt, colors=col[i % 4],levels=[0.5], linewidths=c[i % 4], linestyles='-.')
        if i == 1:
            sR.nYplot=3
            sR.plotContour(axt, levels=[0.5],colors=col[i % 4],  linewidths=c[i % 4], linestyles='-.')
        if i == 2:
            sR.nYplot=3
            sR.plotContour(axt, levels=[0.5], linewidths=c[i % 4], linestyles='-.',colors=col[i % 4])           
        sR.calculateNormVelocity()
        sR.normV[np.where(sR.normV==0)]=np.nan
        contours=d6py.Tools.findContours(sR.grid_x[:,0, 0], sR.grid_z[0,0,:], sR.normV[:,sR.nYplot,:], 0.01)
        if ifile<nF:
            for cont in contours:
                linesC= axs[1, k].plot(smooth(cont[:, 0],5)*100, smooth(cont[:, 1],5)*100, linestyle=ls2[i], color=col[i % 4], linewidth=c[i % 4], alpha=1, label="limit-mod")
  
            
        V = area(sR.findContourPhi(level=0.5)[0])

        lost = (Vini[i]-V)/Vini[i]*100
        print(lost)
        mod = 'velocity'

        if ifile / sR.dConfig['fps'] < 0.5:
            X=np.array([sR.datas[ifile,2],sR.datas[ifile,2]])/sR.scaleLength
            Y=np.array([sR.datas[ifile,4]+sR.Ldoor/2,sR.datas[ifile,4]-sR.Ldoor/2])/sR.scaleLength
            axt.plot([0, 0], Y, 'k', alpha=0.5)

        axt.set_yticklabels([])
        axt.set_xticklabels([])

        if ifile<=nF:
            h = sR.CS.legend_elements(str(sR.dimSim)+"D~-~ \mu= " + toS(sR.dConfig['mu'], 2))[0]
            
            if i ==0:
                l = [r"Sim. free surf. W = "+format(sR.dConfig['box'][1]*100*6/8,'.0f') + ' cm']
            elif i==1:
                l = [r"Sim. free surf. W = "+format(sR.dConfig['box'][1]*100*6/8,'.0f') + ' cm']
                for co in sR.CS.collections:
                    co.set_dashes([(0.0, [4, 2, 1, 2, 1, 2])])
                h = sR.CS.legend_elements(str(sR.dimSim)+"D~-~ \mu= " + toS(sR.dConfig['mu'], 2))[0]
            elif i==2:
                l = [r"Sim. free surf. W = "+format(sR.dConfig['box'][1]*100-4,'.0f') + ' cm']
            h1, l1 = h1+h, l1+l
            if i==0:
                l = [r"Sim. static-flowing W = "+format(sR.dConfig['box'][1]*100*6/8,'.0f') + ' cm']
            elif i==1:
                l = [r"Sim. static-flowing W = "+format(sR.dConfig['box'][1]*100*6/8,'.0f') + ' cm']
            elif i==2:
                l = [r"Sim. static-flowing W = "+format(sR.dConfig['box'][1]*100-4,'.0f') + ' cm']
            h1, l1 = h1+ linesC, l1+l  
            
    sR.plotDoor(ax, alpha=0.5)
    k += 1
    
l_new=[l1[2],l1[0],l1[4],l1[3],l1[1],l1[5]]
h_new=[h1[2],h1[0],h1[4],h1[3],h1[1],h1[5]]
fig.legend(h_new  , l_new , fontsize=leg_fontsize,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=2,bbox_to_anchor=(0.15, 0.01, 0.4, 0.2))

fig.savefig(driveFolder+"/programs/gitLab/dry-granular/doc/article/figures/"+fignames[Nrun]+"_width_var.pdf", dpi=150)
plt.show()
# %%
