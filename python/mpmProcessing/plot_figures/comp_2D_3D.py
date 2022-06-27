#%%
%matplotlib qt5
from default_2_lines_plot import *


#%

Nrun=4
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

# R1_2d, selectedDict = d6py.whereSand6OutFromParms(listNumRun, mu=0.58, delta_mu=0., runNumber=Nrun, dimSim=2, delta_mu_start=0)


# R2_2d, selectedDict = d6py.whereSand6OutFromParms(listNumRun, mu=0.46, delta_mu=0., runNumber=Nrun, dimSim=2, delta_mu_start=0)

R_ref, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=3, delta_mu=0.0,delta_mu_start=0,muRigid=0.23, keyWord='W_6.0_R_1.0')

R2, selectedDict = d6py.whereSand6OutFromParms(listNumRun,  runNumber=Nrun, dimSim=3, delta_mu=0.0,delta_mu_start=0,muRigid=0.23, keyWord='W_1.5_R_1.0')



#%%


# selectedRuns = [R_ref[0],R2[1], R1_2d[-1], R2_2d[-1]] #8
# selectedRuns = [R1_2d[-1], R2_2d[-1]] #8
selectedRuns = [R_ref[1],R2[-1]] #8
for sR in selectedRuns:
    sR.scLength(0.01)
    nF = int(sR.dConfig['nFrames'])
    
    



#%
plt.close('all')

fig, axs = figure_2lines_template()
draw_gravity_2_lines(axs, runExp1.dictE['Slope'],fontsize_g)

k = 0
NsR = len(selectedRuns)
indContrst = 0
ls = ['-.', '-.', '-.', '--']
ls2= [':',':',':',':']
c = [1.2, 1, 0.8, 0.8]
c2 = [0.9, 1. , 0.6]
col = [ 'k', 'darkred', 'yellowgreen','darkblue']
SR = selectedRuns[0:4]
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


for ifile in [3, 8, 30]:
# for ifile in [6, 12, 3]:
    print(ifile)
    ax = axs[0, k]
    h1, l1 = [], [] 

    for sR, i in zip(SR, range(len(SR))):
        axt = axs[1, k]
        sR.loadVTK(int(ifile * sR.dConfig['fps'] / 15))
        sR.plotContour(axt, colors=col[i % 4],levels=[0.5], linewidths=c[i % 4], linestyles=ls[i])
        
        sR.calculateNormVelocity()
        sR.normV[np.where(sR.normV==0)]=np.nan
        if sR.dimSim==3:
            contours=d6py.Tools.findContours(sR.grid_x[:,0, 0], sR.grid_z[0,0,:], sR.normV[:,sR.nYplot,:], 0.01)
        else: 
            contours=d6py.Tools.findContours(sR.grid_x[:,0], sR.grid_y[0, :], sR.normV, 0.01)
        if ifile<=nF:
            for cont in contours:
                linesC= axs[1, k].plot(smooth(cont[:, 0],5)*100, smooth(cont[:, 1],5)*100, linestyle=ls2[i], color=col[i % 4], linewidth=c[i % 4], alpha=1, label="limit-mod")
        # if i==2:

        #     sR.opyfPointCloudColoredScatter(axt, s=2)
        V = area(sR.findContourPhi(level=0.5)[0])

        lost = (Vini[i]-V)/Vini[i]*100

        mod = 'velocity'
        # if i ==1:
        #     if ifile / sR.dConfig['fps'] < 0.5:
        #         X=np.array([sR.datas[ifile,2],sR.datas[ifile,2]])/sR.scaleLength
        #         Y=np.array([sR.datas[ifile,4]+sR.Ldoor/2,sR.datas[ifile,4]-sR.Ldoor/2])/sR.scaleLength
        #         axt.plot([0, 0], Y, 'k', alpha=0.5)

        axt.set_yticklabels([])
        axt.set_xticklabels([])

        if i >1:
            for co in sR.CS.collections:
                co.set_dashes([(0.0, [4, 2, 1, 2, 1, 2])])
        h = sR.CS.legend_elements(str(sR.dimSim)+"D~-~ \mu= " + toS(sR.dConfig['mu'], 2))[0]
        if sR.dimSim == 2:
            l = [r"Num. $\mu_{"+format(sR.dimSim,'1.0f')+"D }=$"+format(sR.dConfig['mu'],'0.2f')]
        if sR.dimSim == 3:
            l = [r"Num. ($W="+format(sR.dConfig['box'][1]*100*6/8,'0.0f')+"~\mathrm{cm})$ " ]
                          
        h1, l1 = h1+ h, l1+l  
    # sR.plotDoor(ax, alpha=0.5)
    k += 1
    
folders = ["B00_H_12_W_1_R1_a",
           "B00_H_12_W_1_R1_b",
           "B00_H_12_W_4_R1_a",
           "B00_H_12_W_4_R1_b",
           "B00_H_13_W_1_R1.1_a",
           "B00_H_13_W_4_R1.1_a"]
datas = []
cs= ['darkblue', 'blue', 'red','darkred', ]
labels=['Exp. $W=1$ cm | A', 'Exp. $W=1$ cm | B', 'Exp. $W=4$ cm | A', 'Exp. $W=4$ cm | B', 'W_1cm_R_1.1_a', 'W_4cm_R_1.1_a']

scale=30
i =0

mainF = "/media/gauthier/DataSSD/TAF/TAF_inria/MPM-data/Collapse_Experiment/Vienna/"

outF='/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/walls/experimental_Hf_xf/'

# import shutil
# for f in folders[:]:
#     shutil.copyfile(mainF+'/'+f+'/final_profile/final_profile.csv',outF+'final_profile_'+f+'.csv')


h2, l2 = [],[]
for f in folders[:]:
    H, data = opyf.read_csv(mainF+'/'+f+'/final_profile/final_profile.csv')
    xf=(data[:-3,5]-200-12*30)/scale
    hf=(-data[:-3,6]+495)/scale
    fileN=outF+'final_profile_'+f+'.csv'
    variables =[['xf [m]',xf/100],['hf [m]',hf/100]]
    opyf.write_csvScalar(fileN,variables)
    line=axt.plot(xf,hf,'-o',c=cs[i], label=labels[i], lw=0.5, ms=0.5,zorder=-1)
    datas.append(data)
    h2, l2 = h2 +line  , l2+  [labels[i]]
    i+=1  
    
  
fig.legend(h1  , l1 , fontsize=leg_fontsize-1,loc=3, framealpha=0.5,edgecolor='w',facecolor='w',ncol=1,bbox_to_anchor=(0.5, 0.72, 0.85, 0.2))
fig.legend(h2  , l2 , fontsize=leg_fontsize-1,loc=3, framealpha=0.,edgecolor='w',facecolor='w',ncol=2,bbox_to_anchor=(0.05, 0.01, 0.4, 0.2))

axt.plot([-12, 0.0,0.0,40],[12,12,0,0],'-',color='seagreen',zorder=-1,lw=1.5)

for ax in axs[1, :]:
    ax.set_ylim(-3,15)
    ax.set_xlim(-15,55)


axs[1, 0].remove()
axs[1,1].remove()
axt.set_position([0.05,0.2,0.9,0.8])    
fig.set_size_inches(3.,1.8)
ax.set_xlim(-11.6,25)
for t in fig.texts:
    t.set_visible(False)
plt.show()
fig.savefig(driveFolder+"/programs/gitLab/dry-granular/doc/article/figures/"+fignames[Nrun]+"_R_1_num_vs_exp_W_1cm_and_4cm.pdf", dpi=1200)
  

# %%
