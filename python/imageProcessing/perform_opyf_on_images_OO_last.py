#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 15:36:55 2019

@author: Gauthier Rousseau
"""
#%%
import os, sys, json
import numpy as np
import cv2
os.chdir('/media/gauthier/Data-Gauthier/programs/gitLab/sand6/python/imageProcessing')
from Tools_collapses import mask_collapses
import matplotlib.pyplot as plt


driveFolder='/media/gauthier/Data-Gauthier/Gauthier/'
# sys.path.append(driveFolder+'TAF/TAF_EPFL/current_work/OPyF-Project/github/opyFlow')
import opyf
sys.path.append(driveFolder+'/TAF_inria/INRIA_current_work/GitLab/dry-granular/data/outputs_experiments_and_mpm/python/imageProcessing')


folder_main='/media/gauthier/Data-Gauthier/Gauthier/TAF/TAF_inria/MPM-data/Collapse_Experiment/Video_src/'
JSONpath=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/granular_collapases_imaging_velocity_fields_and_free_surface_elevation/Granular_Collapses_Experimental_Informations.json'
in_file = open(JSONpath,"r")
dictExp = json.load(in_file) 
lExp=[]
for d in dictExp:
    lExp.append(d)        
lExp=np.sort(lExp) 

#def runOpyfonLandSlides(j):
#j=0
#%%
for j in range(6,7):
    sE=lExp[j] #Selected experiment
    sdictE=dictExp[sE]

    mainOutPutFolder=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs_opyf2/'
    workingFolder=mainOutPutFolder+str('Run_'+format(j,'02.0f')+'_'+sdictE['grainType']+'_Slope='+format(sdictE['Slope'],'.0f')+'deg_H='+format(sdictE['H'],'.3f')+'m_L='+format(sdictE['L'],'.3f'))    
    
    opyf.mkdir2(workingFolder)

    writingExpPath=workingFolder
    writingImgCheckPath=writingExpPath+'/DATA_Img_check'

    opyf.mkdir2(writingExpPath)
    opyf.mkdir2(writingImgCheckPath)
    
    #Create a folder to save outputs (csv files, vtk files, images)
    
    
    vidPath= folder_main+sdictE['videosrc']

#creat the video Analyzer Object
    
    vidA=opyf.videoAnalyzer(vidPath,vlim=[0,30],imageROI=sdictE['ROI'])
    step,shift=8,5
    start_frame=int(sdictE['framedeb'])-100
    Ntot = (sdictE['framFin']+400-start_frame)/shift
    vidA.set_vecTime(starting_frame=start_frame,step=step,shift=shift,Ntot=Ntot)
    
    vidA.set_gridToInterpolateOn(stepGrid=2)

#set Good Features To track Params

    feature_params = dict( maxCorners = 40000,
                               qualityLevel = 0.08,
                               minDistance = 3,
                               blockSize = 10)
    
    vidA.set_goodFeaturesToTrackParams(**feature_params)
    
    lk_params = dict( winSize  = (30, 15),
                      maxLevel = 2,
                      criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 50, 0.03))
    
    vidA.set_opticalFlowParams(**lk_params)
    
    filters_params = dict(  RadiusF=20.*sdictE['scale'],
                          minNperRadius=2,
                          maxDevInRadius=np.inf,
                          wayBackGoodFlag=10.,
                          CLAHE=True)
    
    vidA.set_filtersParams(**filters_params)
    
#    vidA. set_interpolationParams(Radius=40, Sharpness=8, kernel='Gaussian')

    
    OR2=sdictE['OR2']
    BW,vecXE,HE=mask_collapses(opyf.Tools.convertToGrayScale(vidA.cropFrameInit),OR2,0.8,sdictE['scale'])
    OR2[1]=OR2[1]-np.mean(HE[-int(len(HE)/2):-30])/sdictE['scale']
    

    vidA.scaleData(metersPerPx=sdictE['scale'], framesPerSecond=sdictE['fps'], origin=OR2)
    vidA.set_vlim([0,1])
    interp_params= dict (Radius=10*sdictE['scale'], # it is not necessary to perform unterpolation on a high radius since we have a high number of values
                 Sharpness=8,
                 kernel='Gaussian')

    vidA.set_interpolationParams(**interp_params)
    vidA.reset()
    HEfinal=[]
    Xaccu = np.empty((0, 2))
    Vaccu = np.empty((0, 2))
    kAccu,stepAccu=0,4
    for pr,i in zip( vidA.prev, vidA.vec):
        kAccu+=1
        vidA.stepGoodFeaturesToTrackandOpticalFlow(pr, i)    

        vidA.mask,vecXE,HE=mask_collapses(opyf.Tools.convertToGrayScale(vidA.vis),OR2,sdictE['maskTh'],sdictE['scale'])

        BWresized = cv2.resize(vidA.mask,(vidA.Lgrid,vidA.Hgrid))  
        vidA.gridMask=np.ones(BWresized.shape)
        vidA.gridMask[np.where(BWresized==0)]=np.nan
        Xaccu = np.append(Xaccu, vidA.X, axis=0)
        Vaccu = np.append(Vaccu, vidA.V, axis=0)
        if pr==True and kAccu==stepAccu:      
            HEfinal.append(HE)
            vidA.Xdata.append(Xaccu)
            vidA.Vdata.append(Vaccu)
            vidA.interpolateOnGrid(Xaccu,Vaccu)
            vidA.UxTot.append(np.reshape(vidA.interpolatedVelocities[:,0],(vidA.Hgrid,vidA.Lgrid)))
            vidA.UyTot.append(np.reshape(vidA.interpolatedVelocities[:,1],(vidA.Hgrid,vidA.Lgrid)))                            
            Field=opyf.Render.setField(vidA.Ux,vidA.Uy,'norm')

            vidA.opyfDisp.plotField(Field,vis=vidA.vis)

            vidA.opyfDisp.ax.plot(vecXE,HE,'--k',linewidth=0.5)
            vidA.opyfDisp.fig.savefig(writingImgCheckPath+'/Field_'+format(i,'04.0f')+'_to_'+format(i+vidA.paramVecTime['step'],'04.0f')+'.png')
#           

#            plt.pause(0.01)
            Xaccu = np.empty((0, 2))
            Vaccu = np.empty((0, 2))
            kAccu=0

    
    vidA.Time=(vidA.Time[0:-1:2]+vidA.Time[1::2])/2
    HEfinal=np.array(HEfinal)    
    
    filename='Run_'+format(j,'02.0f')+'_unstructured_good_features_to_track_velocities'
    vidA.writeGoodFeaturesPositionsAndDisplacements(outFolder=writingExpPath,filename=filename)
#    
    filename='Run_'+format(j,'02.0f')+'_rectilinear_velocity_fields'
    vidA.fieldResults='time-serie'
    vidA.writeVelocityField(outFolder=writingExpPath,filename=filename)
#       
#    #    Export flow depth
    filename=writingExpPath+'/Run_'+format(j,'02.0f')+'_depths.hdf5'

    opyf.hdf5_Write(filename,[['T [s]',vidA.Time],['X [m]',vecXE]],[['Depth [m]',HEfinal]])
    opyf.hdf5_Read(filename)


# %%
