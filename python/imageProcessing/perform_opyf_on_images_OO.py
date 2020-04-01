#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 15:36:55 2019

@author: Gauthier Rousseau
"""

import os, sys, json
import numpy as np
import cv2
from Tools_collapses import mask_collapses


driveFolder='./'


import opyf

#If not in the repository, you may download it on https://github.com/groussea/opyflow


vidPath='/media/gauthier/Samsung_T5/MPM_data/Collapse_experiment/Video_src'

driveFolder='/media/gauthier/Data-Gauthier/Gauthier/'
folder_main=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Video_src' 
os.chdir(folder_main)
JSONpath=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Video_src/dictExp.json'
in_file = open(JSONpath,"r")
dictExp = json.load(in_file) 
lExp=[]
for d in dictExp:
    lExp.append(d)        
lExp=np.sort(lExp) 

#def runOpyfonLandSlides(j):
#j=0
#%%
for j in range(1,8):
    sE=lExp[j] #Selected exeperiment
    sdictE=dictExp[sE]
    mainOutPutFolder=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs_opyf2/'
    workingFolder=mainOutPutFolder+str('Run_'+format(j,'02.0f')+'_'+sdictE['grainType']+'_Slope='+format(sdictE['Slope'],'.0f')+'deg_H='+format(sdictE['H'],'.3f')+'m_L='+format(sdictE['L'],'.3f'))    
    
    opyf.mkdir2(workingFolder)
    writingEDProfPath=workingFolder+'/DATA_Exp_Depth'
    writingTracksPath=workingFolder+'/DATA_Exp_Tracks'
    writingImgCheckPath=workingFolder+'/DATA_Img_check'
    opyf.mkdir2(writingEDProfPath)
    opyf.mkdir2(writingTracksPath)
    opyf.mkdir2(writingImgCheckPath)
    
    #Create a folder to save outputs (csv files, vtk files, images)
    
    
    vidPath= folder_main+sdictE['videosrc']

#creat the video Analyzer Object
    
    vidAnalyzer=opyf.videoAnalyzer(vidPath,vlim=[0,30],imageROI=sdictE['ROI'])
    step,shift,Ntot=8,10,int(sdictE['nFrames'])+5
    vidAnalyzer.set_vecTime(starting_frame=int(sdictE['framedeb'])-100,step=step,shift=shift,Ntot=Ntot)
    
    vidAnalyzer.set_gridToInterpolateOn(stepGrid=4)

#set Good Features To track Params

    feature_params = dict( maxCorners = 40000,
                               qualityLevel = 0.08,
                               minDistance = 3,
                               blockSize = 10)
    
    vidAnalyzer.set_goodFeaturesToTrackParams(**feature_params)
    
    lk_params = dict( winSize  = (30, 15),
                      maxLevel = 2,
                      criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 50, 0.03))
    
    vidAnalyzer.set_opticalFlowParams(**lk_params)
    
    filters_params = dict(  RadiusF=20.,
                          minNperRadius=2,
                          maxDevInRadius=np.inf,
                          wayBackGoodFlag=10.)
    
    vidAnalyzer.set_filtersParams(**filters_params)
    

    
    OR2=sdictE['OR2']
    BW,vecXE,HE=mask_collapses(opyf.Tools.convertToGrayScale(vidAnalyzer.cropFrameInit),OR2,sdictE)
    OR2[1]=OR2[1]-np.mean(HE[-int(len(HE)/2):-30])/sdictE['scale']
    

    vidAnalyzer.scaleData(metersPerPx=sdictE['scale'], framesPerSecond=sdictE['fps'], origin=OR2)
    vidAnalyzer.set_vlim([0,1])
    interp_params= dict (Radius=10*sdictE['scale'], # it is not necessary to perform interpolation on a high radius since we have a high number of values
                 Sharpness=8*sdictE['scale'],
                 kernel='Gaussian')

    vidAnalyzer.set_interpolationParams(**interp_params)
    vidAnalyzer.reset()
    HEfinal=[]
    for pr,i in zip( vidAnalyzer.prev, vidAnalyzer.vec):
        vidAnalyzer.stepGoodFeaturesToTrackandOpticalFlow(pr,i)  
        vidAnalyzer.mask,vecXE,HE=mask_collapses(opyf.Tools.convertToGrayScale(vidAnalyzer.vis),OR2,sdictE)
        BWresized = cv2.resize(vidAnalyzer.mask,(vidAnalyzer.Lgrid,vidAnalyzer.Hgrid))  
        vidAnalyzer.gridMask=np.ones(BWresized.shape)
        vidAnalyzer.gridMask[np.where(BWresized==0)]=np.nan

        if pr==True:      
            HEfinal.append(HE)
            vidAnalyzer.Xdata.append( vidAnalyzer.X)
            vidAnalyzer.Vdata.append(vidAnalyzer.V)
            vidAnalyzer.interpolateOnGrid(vidAnalyzer.X,vidAnalyzer.V)
            vidAnalyzer.UxTot.append(np.reshape(vidAnalyzer.interpolatedVelocities[:,0],(vidAnalyzer.Hgrid,vidAnalyzer.Lgrid)))
            vidAnalyzer.UyTot.append(np.reshape(vidAnalyzer.interpolatedVelocities[:,1],(vidAnalyzer.Hgrid,vidAnalyzer.Lgrid)))                            
            Field=opyf.Render.setField(vidAnalyzer.Ux,vidAnalyzer.Uy,'norme')
#            if display=='field':
            vidAnalyzer.opyfDisp.plotField(Field,vis=vidAnalyzer.vis)
            vidAnalyzer.opyfDisp.fig.savefig(writingImgCheckPath+'/Field_'+format(i,'04.0f')+'_to_'+format(i+vidAnalyzer.paramVecTime['step'],'04.0f')+'.png')
#
#                    if saveImgPath is not None:
#                        self.opyfDisp.fig.savefig(saveImgPath+'/'+display+'_'+format(i,'04.0f')+'_to_'+format(i+self.paramVecTime['step'],'04.0f')+'.'+imgFormat)
#                    plt.pause(0.01)
            
#            plt.pause(0.01)
    
    
    filename='Run_'+format(j,'02.0f')+'_unstructured_good_features_to_track_velocities'
    vidAnalyzer.writeGoodFeaturesPositionsAndDisplacements(outFolder=writingTracksPath,filename=filename)
    
    filename='Run_'+format(j,'02.0f')+'_rectilinear_velocity_fields'
    vidAnalyzer.writeVelocityField(outFolder=writingTracksPath,filename=filename)
       
    #    Export flow depth
    filename=writingEDProfPath+'/Run_'+format(j,'02.0f')+'_depths.hdf5'
    HEfinal=np.array(HEfinal)
    opyf.hdf5_Write(filename,[['T [s]',vidAnalyzer.Time],['X [m]',vecXE]],[['Depth [m]',HEfinal]])
    opyf.hdf5_Read(filename)
#    vidAnalyzer.extractGoodFeaturesPositionsAndDisplacements(display='quiver',displayColor=True,nomalize=False,s=10,nvec=5000)

#%%    


#%%


