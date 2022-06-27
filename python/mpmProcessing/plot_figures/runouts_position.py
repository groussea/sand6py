
#%%
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


#%
%matplotlib qt

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth  
plt.close('all')
fig, ax, ax12= generate_fig()
fig2, ax2, ax22= generate_fig()

mainExpFolder = driveFolder + \
        '/TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/'
cmap=plt.get_cmap('inferno')
ic=0
npl=1
hth,vth=0.007, 0.001 # Threshold on runout height in (m) and on final time en m/s
tfsM,tfsE=[],[]
outf='/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/transient_runouts_B_runs/'
for Nrun in range(4,4+npl):
    run=fignames[Nrun]
    R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0, runNumber=Nrun, dimSim=3, delta_mu_start=0.0, muRigid=0.23, mu=0.44,keyWord='W_8')
    xfMss, VxfMss, tfs, maxVss, times, dt= loadRunOutMod(R1[0],hth=hth,vth=vth)
    tfsM.append(tfs[0])
    VxfMssP=np.insert(VxfMss[0][:-1:],0,0)
    timesP=np.insert(times[:-2:]+dt/2,0,0)
    plot=ax.plot(timesP, VxfMssP,'-',c=cmap(ic/npl+0.05),lw=0.9, label=run)

    runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
    xfEs, timesE, maxVsE, VxfEs, tf_exp,timesE, dtE = loadRunOutExp(runExp1,hth=hth,vth=vth)
    tfsE.append(tf_exp)
    VxfMssEP=np.insert(smooth(VxfEs[10:-20:],25),0,0)
    timesEP=np.insert(timesE[10:-21:]+dtE/2,0,0)
    plotE=ax.plot(timesEP, VxfMssEP,':',c=cmap(ic/npl+0.05),lw=0.9)
    ax2.plot(times[:-2], xfMss[0][:-2],'-',c=cmap(ic/npl+0.05),lw=0.9, label=' '+run)
    ax2.plot(timesE[8:-10], xfEs[8:-10],':',c=cmap(ic/npl+0.05),lw=0.9)
    opyf.mkdir2(outf+run)
    filename1=outf+run+'/'+run+'_exp_runout_pos_Hth='+format(hth*1000,'1.0f')+'_mm.csv'
    filename2=outf+run+'/'+run+'_mod_runout_pos_Hth='+format(hth*1000,'1.0f')+'_mm.csv'
    filename3=outf+run+'/'+run+'exp_runout_vel_Hth='+format(hth*1000,'1.0f')+'_mm.csv'
    filename4=outf+run+'/'+run+'mod_runout_vel_Hth='+format(hth*1000,'1.0f')+'_mm.csv'
    
    variables1 = [['time [s]',timesE[8:-10]],['exp_runout_pos [m]',xfEs[8:-10]]]
    variables2 = [['time [s]',times[:-2]],['mod_runout_pos [m]', xfMss[0][:-2]]]
    variables3 = [['time [s]',timesEP],['exp_runout_vel [m/s]', VxfMssEP]]
    variables4 = [['time [s]',timesP],['mod_runout_vel [m/s]',VxfMssP]]
    
    opyf.write_csvScalar(filename1, variables1)
    opyf.write_csvScalar(filename2, variables2)
    opyf.write_csvScalar(filename3, variables3)
    opyf.write_csvScalar(filename4, variables4)
    ic+=1


ax.legend(fontsize=leg_fontsize,loc=5, framealpha=1,edgecolor='w',facecolor='w',ncol=1,bbox_to_anchor=(0.6, 0.7, 0.4, 0.2))
ax2.legend(fontsize=leg_fontsize,loc=5, framealpha=1,edgecolor='w',facecolor='w',ncol=1,bbox_to_anchor=(0.6, 0.1, 0.4, 0.2))    
ax.set_ylabel('d$L_t$/d$t$ [m/s]')
ax2.set_ylabel('$L_t-L_0$ [m]')

fig.text(0.03,0.94,'(a)',fontsize=9)
fig2.text(0.03,0.94,'(b)',fontsize=9)

ax.set_xlim([0, 1.2])
fig.show()
fig2.show()

#%% save runouts data in csv





#%%
os.chdir('/media/gauthier/DataSSD/programs/gitLab/dry-granular/doc/article/figures/front_vel/')
fig.savefig('B_fronts_vel.pdf')
fig2.savefig('B_fronts_pos.pdf')
# ax.plot(times[:-2], xfMss[0][:-2],'+:',c='g',lw=0.9, label='Num. front pos. '+run)
# # plt.plot(times[:-2], xfMss[1][:-2],'1--',c='y',lw=1.2, label='Num. front position '+run+' visc')
# ax.plot(timesE[8:-10], xfEs[8:-10],'--',c='k',lw=1.4, label='Exp. front pos. '+run)

# %% save in csv final runouts and heigths

plt.close('all')
import csv 

fig, ax =plt.subplots(1,1)

H0s,L0s=[],[] 
xfs,Hfs=[],[]
xfsE,HfsE=[],[]
tfsE, tfsM=[],[]
vth=0.001

filename = "/media/gauthier/DataSSD/programs/gitLab/dry-granular/data/final_runouts_heights_tf_mod_vs_exp/final_runouts_heights_tf_mod_vs_exp_all_BG_runs_Hth=1mm_Vth="+format(vth,'.3f')+"ms-1.csv"

for Nrun in range(0,9):
    run=fignames[Nrun]
    if Nrun<4:
        muR, mu = 0.3, 0.75
    else:
        muR, mu = 0.23, 0.44
    R1, selectedDict = d6py.whereSand6OutFromParms(listNumRun, delta_mu=0, runNumber=Nrun, dimSim=3, delta_mu_start=0.0, muRigid=muR, mu=mu,keyWord='W_8')
    xfMss, VxfMss, tfs, maxVss, times, dt= loadRunOutMod(R1[0],hth=hth,vth=vth)
    sR=R1[0]
    nF = int(sR.dConfig['nFrames'])
    xfM, tf=extractRunOut_time(sR,nF,th=hth)
    Hfs.append(np.max(sR.pointsp[:,2]))
    xfs.append(xfM)
    runExp1 = d6py.ExperimentalRun(Nrun, mainExpFolder, loadField=True)
    H0s.append(runExp1.dictE['H'])
    L0s.append(runExp1.dictE['L'])
    xfEs, timesE, maxVsE, VxfEs, tf_exp,timesE, dtE = loadRunOutExp(runExp1,hth=hth,vth=vth)
    HfsE.append(np.max(runExp1.expD[-1][30:]))
    xfsE.append(np.nanmax(xfEs))
    tfsE.append(tf_exp)
    tfsM.append(tfs[0])
    ax.plot(runExp1.expD[-1])
    
print(Hfs)
print(HfsE)

fig.show()
variables = [['Run',fignames],['H0_[m]',H0s],['Hf_exp_[m]',HfsE],['Hf_mod_[m]',Hfs],['L0_[m]',L0s],['xf_exp_[m]',xfsE],['xf_mod_[m]',xfs],['tf_exp_[m]',tfsE],['tf_mod_[s]',tfsM]]

opyf.write_csvScalar(filename, variables)
