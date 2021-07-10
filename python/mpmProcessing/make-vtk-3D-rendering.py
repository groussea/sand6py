
# %%
#!/usr/bin/env python3
import sys
import os
import opyf  # from opyflow library some rendering function may be employed
import numpy as np
sys.path.append(
    '/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python/imageProcessing')
from Tools_collapses import mask_collapses
import matplotlib.pyplot as plt
from d6py.Tools import *
import d6py
# intialize exteral packages

import matplotlib

maind6OutFolder = '/media/gauthier/Samsung_T5/sand6_sorties/sand6_out'
paths, folders, listDictConf, listNumRun = d6py.findOutSand6Paths(maind6OutFolder, 4)

# d6py.whereSand6OutFromFolder(listNumRun,maind6OutFolder )
# N = 7
# Reference
# Ref, selectedDict = d6py.whereSand6OutFromParms(listNumRun,fps=150, runNumber=N, dimSim=3)
# Ref, selectedDict = d6py.whereSand6OutFromParms(listNumRun,res=[122.0, 6.0, 30.0],substeps=2,muRigid=0.18,mu=0.38, delta_mu=0, runNumber=N, dimSim=3)

# path='/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/rotating_drum/rot_drum_mu=0.75_muRigid=2_delta_mu=0.000_substeps_10_I0_start=0.0040_delta_mu_start=0.0000'
path = '/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/rotating_drum/rot_drum_mu=0.75_muRigid=2_delta_mu=0.000_substeps_25_I0_start=0.0000_delta_mu_start=0.00002'
path='/media/gauthier/Samsung_T5/sand6_sorties/sand6_out/Run_04_3D_Door_mu=0.48_muRigid=0.05_W_8.0cm_glass-beads-0.47mm_Slope=0deg_delta_mu=0.00_substeps_4_fracH=0.8_I0_start=0.0040_delta_mu_start=0.00cohesion_render/'
sR = d6py.NumericalRun(path)
nF = int(sR.dConfig['nFrames'])


from pyevtk.hl import pointsToVTK, gridToVTK
outputPath='/media/gauthier/Samsung_T5/sand6_sorties/for_rendering/collapse_with_co/'


#%%
k = 0
sR.loadVTK(k)
(l,c)=sR.vel.shape
values = np.zeros((l, c + 1))



#%%
for k in range(0,nF):
    sR.loadVTK(k)
    path = outputPath + '/points_' + format(k, '04.0f')
    values[:, 0:3] = sR.vel
    values[:, 3] = sR.vol
    
    # Interp=npOutputValues=opyf.Interpolate.npInterpolateVTK3D(sR.pointsp,values,Tpoints,ParametreInterpolatorVTK=ParametreInterpolatorVTK)

    x=np.ascontiguousarray(sR.pointsp[:,0],dtype=np.float32)
    y=np.ascontiguousarray(sR.pointsp[:,1],dtype=np.float32)
    z=np.ascontiguousarray(sR.pointsp[:,2],dtype=np.float32)
    vol = np.ascontiguousarray(sR.vol, dtype=np.float32)
    velx = np.ascontiguousarray(sR.vel[:, 0], dtype=np.float32)
    vely = np.ascontiguousarray(sR.vel[:, 1], dtype=np.float32)
    velz = np.ascontiguousarray(sR.vel[:, 2], dtype=np.float32)
    inertia = np.ascontiguousarray(sR.inertia, dtype=np.float32)
    pointsToVTK(path, x, y, z, data = {"vol" : sR.vol,"velx" : velx,"vely" : vely,"velz" : velz,"inertia" :inertia})






# ParametreInterpolatorVTK = {'kernel': 'Gaussian',
#                             'Radius': 0.0035,
#                             'Sharpness': 2.}
    #del Tpoints
    # Ux=np.ascontiguousarray(Interp[:,0],dtype=np.float32)
    # Uy=np.ascontiguousarray(Interp[:,1],dtype=np.float32)
    # Uz = np.ascontiguousarray(Interp[:, 2], dtype=np.float32)
    # vol = np.ascontiguousarray(Interp[:, 2], dtype=np.float32)
    
    # Vx=opyf.Interpolate.npTargetPoints2Grid3D(Interp[:,0],resP[0],resP[1],resP[2])
    # Vy=opyf.Interpolate.npTargetPoints2Grid3D(Interp[:,1],resP[0],resP[1],resP[2])
    # Vz=opyf.Interpolate.npTargetPoints2Grid3D(Interp[:,2],resP[0],resP[1],resP[2])
    # vol=opyf.Interpolate.npTargetPoints2Grid3D(Interp[:,3],resP[0],resP[1],resP[2])

    # dx, dy, dz=vecX[1]-vecX[0],vecY[1]-vecY[0],vecZ[1]-vecZ[0]

    # vecX_cell=np.linspace(vecX[0]-dx/2,vecX[-1]+dx/2,int(resP[0]+1))
    # vecY_cell=np.linspace(vecY[0]-dy/2,vecY[-1]+dy/2,int(resP[1]+1))
    # vecZ_cell=np.linspace(vecZ[0]-dz/2,vecZ[-1]+dz/2,int(resP[2]+1))

    # grid_xC, grid_yC, grid_zC= np.meshgrid(vecX_cell,vecY_cell,vecZ_cell)

    # path = outputPath + '/structured_150_' + format(k, '04.0f')
    # gridToVTK(path, grid_xC, grid_yC,  grid_zC,cellData = {"Ux" : Vx, "Uy" : Vy, "Uz" : Vz, "vol" : vol})
