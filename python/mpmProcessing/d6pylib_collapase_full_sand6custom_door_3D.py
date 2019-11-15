#!/usr/bin/env python
"""
Created on Mon Oct 16 14:35:06 2017

This python script tranform the Particles-###.vtk files in to csv files
+read the config file to store simulation parameters in a dict
+perform the scaling to obtain real length an velocities
+interpolate the data to obtain the velocity field
+plotting functions
4 Functions can be imported extractInfoVTK, readConfigFile, csv_write, csv_read,
interpolate2D, d6imshow


extractInfoVTK(vtkfilename) return the velocities, the volumes and points of each MPM point
readconfigFile(configFile) return a disctionary with 18 firs lines of the config file 
.... continuing

Usage :
d6pylib.py [writing path] [d6vtk path(by default 'out')]

Warning : only used for 2D simulation (modification are needed for 3D)
If no argments: the csv files are saved in the out folder considering the python file is in the build folder


@author: Gauthier
"""
import cv2
from shutil import copyfile
import os, shutil
import numpy as np  
import csv
import matplotlib.pyplot as plt
import json
import sys
import subprocess
from matplotlib_scalebar.scalebar import ScaleBar

driveFolder='/scratch/garousse/'
driveFolder='/media/gauthier/Gauthier_Backup/'
sys.path.append(driveFolder+'TAF/TAF_EPFL/OPyF-Project/github/opyFlow')
import opyf

sys.path.append(driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/python-essentials/Essentials OpenCV')
from TrackandInterpolate import *
from vtk.util import numpy_support as VN
#import trackandinterpolate

    
    #the path where all the videos are

vidPath=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Video_src' 

    # the genarated dictionnary-json path
JSONpath=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Video_src/dictExp.json'

    # the d6 soft path
d6Path=driveFolder+'TAF/TAF_inria/Sand6/epfl_lhe_2d_and_3d/build_fast'
d6Path='/media/gauthier/Data-Gauthier/programs/gitLab/sand6/build'
#d6Path=driveFolder+'TAF/TAF_inria/GitLab/sand6cohesive/build_julien'
#d6Path='/home/gauthier/programs/epfl_lhe/build2d'
d6OutFolder='out'
out_Opyf=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs_opyf/'
listOutOpyf=['Run_00_gravels-2.7mm_Slope=0deg_H=0.114m_L=0.187',
'Run_01_gravels-2.7mm_Slope=5deg_H=0.124m_L=0.198',
'Run_02_gravels-2.7mm_Slope=10deg_H=0.115m_L=0.279',
'Run_03_gravels-2.7mm_Slope=15deg_H=0.117m_L=0.268',
'Run_04_glass-beads-0.47mm_Slope=0deg_H=0.119m_L=0.223',
'Run_05_glass-beads-0.47mm_Slope=5deg_H=0.111m_L=0.223',
'Run_06_glass-beads-0.47mm_Slope=10deg_H=0.111m_L=0.224',
'Run_07_glass-beads-0.47mm_Slope=15deg_H=0.111m_L=0.224']

os.chdir(d6Path)
#%%

    #add d6py module pythonpath
sys.path.append(d6Path+'/../python')

import d6py
#%
# Load the contents from the file, which creates a new dictionary with all experimental infos
in_file = open(JSONpath,"r")
dictExp = json.load(in_file)   


#generate the list of avi files
lExp=[]
for d in dictExp:
    lExp.append(d)        
lExp=np.sort(lExp) 

import time
t=time.time()

door='with'

for j in range(7,8 ): 
#    plt.close('all')
#for j in range(0,8):
    sE=lExp[j] #Selected exeperiment
    sdictE=dictExp[sE]
    fracH=0.9
    
    if sdictE['camType']=='Phantom':
        sdictE['L']=sdictE['L']+0.1
        sdictE['Ltot']=sdictE['Ltot']+0.1
    #We Consider a 1 meter simu
    L=sdictE['L']    
    Lmod=sdictE['Ltot']
    nFrames=sdictE['nFrames']
    
    if (sdictE['camType']=='BW') & (sdictE['Slope']==15.):
        Lmod=sdictE['Ltot']+0.5
    

    muRigid=0.18
    
    Hmod=(sdictE['H']+0.005)/fracH
    delta_mu=0.
    if door=='with':
        ts=0
    else:
        ts=20
    substeps=4
    sdictE['delta_mu']=delta_mu

    mu=sdictE['mu']
    prop='low'
    if prop=='low':
        mu=sdictE['mu']-0.03
    prop=''
    if prop=='inscrit':
        mu=mu/(1+1/3*mu**2)

        
    if door=='with':
        runName=str('Run_'+format(j,'02.0f')+'_3D_Door_muRigid='+str(muRigid)+'_'+sdictE['grainType']+'_Slope='+format(sdictE['Slope'],'.0f')+'deg_H='+format(sdictE['H'],'.3f')+'m_L='+format(sdictE['L'],'.3f')+'_delta_mu='+format(delta_mu,'.3f')+'_substeps_'+str(substeps)+'_fracH='+str(fracH)+prop+'mu_s='+str(mu))    
    else:
        runName=str('Run_'+format(j,'02.0f')+'_3D_no_Door_start_at_0.13_s_'+sdictE['grainType']+'_Slope='+format(sdictE['Slope'],'.0f')+'deg_H='+format(sdictE['H'],'.3f')+'m_L='+format(sdictE['L'],'.3f')+'_delta_mu='+format(delta_mu,'.3f')+'_substeps_'+str(substeps)+'_fracH='+str(fracH)+prop+'mu_s='+str(mu))      
  
    
     
    d6OutFolder=driveFolder+'TAF/TAF_inria/MPM-data/Collapse_Experiment/Sand6Out/outputs/Tests/'+runName
    opyf.mkdir2(d6OutFolder) 
    
    newConfigFile=d6OutFolder+'/collapse.3d_'+runName+'m.conf'
    configFilein=d6Path+'/../scenes/collapse.3d.LHE.Door.conf'

    d6py.modifyConfigFile(configFilein,newConfigFile,'box',[Lmod, 0.06,Hmod])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'fps',[150])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'gravity',[+9.81*np.sin(sdictE['Slope']*np.pi/180), 0,-9.81*np.cos(sdictE['Slope']*np.pi/180)])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'nFrames',[nFrames+ts])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'randomize',[0])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'substeps',[substeps])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'grainDiameter',[sdictE['grainDiameter']])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'mu',[mu])
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'muRigid',[muRigid]) #mu Door
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'delta_mu',[delta_mu])
    if door=='with':
        if j<=3:
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.06 veldoor:0.8 ts:'+format(ts,'1.0f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L))
        else:
            d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.13 veldoor:0.7 ts:'+format(ts,'1.0f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L))
    else:
        d6py.modifyConfigFile(newConfigFile,newConfigFile,'scenario','collapselhedoor taudoor:0.0001 veldoor:100 ts:'+format(ts,'1.0f')+' frac_h:'+format(fracH,'1.1f')+' column_length:'+str(L))

    
    TypicalLength=0.005
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'res',[Lmod//TypicalLength,0.06//TypicalLength,np.round(Hmod/TypicalLength/10)*10]) #pour avoir un réolution divisible par 10 selon Y
    d6py.modifyConfigFile(newConfigFile,newConfigFile,'I0',[0.279]) 
      
    #load the final config file dictionnary    
        
    dConfigmod=d6py.readConfigFile(newConfigFile)
    #% 

    ################# Execute Different sytem commands ###############     
     
    #proc =subprocess.Popen('rm -rf  ' + d6OutFolder +'/*' , shell=True,stdout=subprocess.PIPE)
    #(out, err) = proc.communicate()    
    cmd = str('./apps/d6 ' +d6OutFolder+ ' -i '+ newConfigFile)
    
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    print("program output:", out)
    
    
    if out=='':
        print('No sand6 simulation were executed - run aborted')
        sys.exit()  
    else:
        print("program output:", out) 
    
    cmd = str('./apps/d62vtk '+d6OutFolder+ ' -a -p')
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
     
    
    if out=='':
        print('No vtk files generated - run aborted')
        sys.exit() 
    else:
        print("program output:", out)
#    shutil.copyfile(d6Path+'/out/door.txt',d6OutFolder+'/door.txt')
    #%%
    ################# Write the Different subfolders ###############     
        
     
    
#    configFileout=d6OutFolder+'/config' 
#    copyfile(configFileout, str(d6OutFolder+'/configout_'+runName))
#    
#    dConfig=d6py.readConfigFile(configFileout)
#        
#       
#        
#    #    f=open(workingFolder+'/info.txt','w')
#    #    writer = csv.writer(f,delimiter=' ')       
#    #    writer.writerow('The simulations begin at the frame ' + str(frameref-framedeb)) 
#    #    f.close()     
#        
#    #Create the grid to interpolate MPM datas
#    res=dConfig['gridScale']/dConfig['nSamples']
#    Nx=int(float(dConfig['box'][0])/res)
#    Ny=int(float(dConfig['box'][1])/res)
#    grid_y, grid_x = np.mgrid[0:int(Ny), 0:int(Nx)]*res
#    
#    
#    #%
#    #    ifile=0
#    
#    writingTracksPath=out_Opyf+listOutOpyf[j]+'/DATA_Exp'
#    filename=writingTracksPath+'/Run_'+format(j,'02.0f')+'_rectilinear_velocity_fields.hdf5'   
#    coordinates,variables=opyf.hdf5_Read(filename)
#    writingEDProfPath=out_Opyf+listOutOpyf[j]+'/DATA_Exp'
#    filename=writingEDProfPath+'/Run_'+format(j,'02.0f')+'_depths.hdf5'    
#    coordinatesDP,variablesDP=opyf.hdf5_Read(filename)
#    
#    
#    
#    UxExp = variables[0][1]
#    UyExp = variables[1][1]
#    UxExp[np.where(UxExp==0)]=np.nan
#    UyExp[np.where(UxExp==0)]=np.nan
#    
#    #    Time=Time-Time[0]
#    Ntot=nFrames+ts
#    shift= 10
#    framedeb=int(sdictE['framedeb'])-100
#    deltaN=Ntot*shift
#    vec1=np.arange(int(framedeb),int(framedeb+deltaN),shift)  
#    cap = cv2.VideoCapture(vidPath+sdictE['videosrc']) 
#    Ratio=(dConfig['columnLength']/(dConfig['box'][0]))
#    Ory=float(dConfig['box'][1])*0.1
#    Orx=float(dConfig['box'][0])*(Ratio)   
#    header, datas=opyf.read_csv(d6OutFolder+"/door.txt")
#    posDoor=datas[:,2:4]-np.array([Orx,Ory])
#    velDoor=datas[:,4:6]
#    res=dConfig['gridScale']/dConfig['nSamples']
#    Nx=int(float(dConfig['box'][0])/res)
#    Ny=int(float(dConfig['box'][1])/res)
#    grid_y, grid_x = np.mgrid[0:int(Ny), 0:int(Nx)]*res
#    
#    
#    grid_xp=grid_x-Orx
#    grid_yp=grid_y-Ory
#    ImgCheckDepth=d6OutFolder+'/ImgCheck_Depth'
#    opyf.mkdir2(ImgCheckDepth)
#    writingMPMDPath=d6OutFolder+'/DATA_Depth'
#    opyf.mkdir2(writingMPMDPath)
#    
#    plt.close('all')   
#    fig,ax=plt.subplots()  
#    shiftMPM=25
#    ifile=0
#    ROI=sdictE['ROI']
#    OR2=sdictE['OR2']
#    Ratioxr=np.float(OR2[0])/ROI[2]
#    Ratioyr=np.float(ROI[3]-OR2[1])/ROI[3]                
#    extentr=np.array([-Ratioxr*ROI[2],(1-Ratioxr)*ROI[2],-Ratioyr*ROI[3],(1-Ratioyr)*ROI[3]])*sdictE['scale']
#    writingcsvPath= d6OutFolder+'/DATA_MPM_particles'
#    opyf.mkdir2(writingcsvPath)
#    d6Depth=[]
#    while ifile < np.floor(nFrames):
#    #while ifile < 1:
#        VTKfilename= d6OutFolder+'/vtk/particles-'+ str(ifile) +'.vtk'        
#        npvelArr,npvolArr,nppointsArr=d6py.extractInfoVTK(VTKfilename)
#        vel=npvelArr*dConfig['velScale']
#        vol=npvolArr*(dConfig['gridScale'])**3
#        points=nppointsArr*dConfig['gridScale']
#        
#        csvFilePath=writingcsvPath+'/particles'+format(ifile,'03.0f')+'.csv'
#     
#        d6py.csv_write(points,vel,vol, csvFilePath)
#        names,values=opyf.read_csv(csvFilePath)
#        points=values[:,1:4]
#        vel=values[:,4:7]
#        vol=values[:,7]
#    
#        pointsp=np.zeros_like(points)
#        pointsp[:,0]=points[:,0]-Orx
#        pointsp[:,1]=points[:,1]-Ory                
#        ####translate the origine on the points and the grid (comment if isnot dConfig['columnLength'])
#        #Write the MPM interpolated filed values in csv files         
#    
#    
#        #% evaluate and save the depth and create a mask
#        hstep=dConfig['box'][0]/dConfig['res'][0]
#        resX=hstep/dConfig['nSamples']*2
#        vecX,H=d6py.depthProfile(grid_xp,pointsp,resX)
#     
#        #Save depth profile
#    
#        csvProf=writingMPMDPath+'/DP_'+format(ifile,'04.0f')+'.csv'
#    
#        #Write the the scaled values in csv files         
#        f=open(csvProf,'w')
#        writer = writer = csv.writer(f, delimiter=';')       
#        writer.writerow(vecX) 
#        writer.writerow(H) 
#        f.close()  
#        
#        h, w=grid_yp.shape
#        mask=np.ones_like(grid_xp)*255
#        
#        for i in range(0,w):    
#            indH=np.argmin(np.absolute(grid_xp[0,i]-vecX))
#            heni=H[indH]
#            indM=np.where(grid_yp[:,i]>heni)
#            mask[indM,i]=0
#        
#        #%
#        interp_params= dict (Radius=dConfig['gridScale']*2/3, # it is not necessary to perform unterpolation on a high radius since we have a high number of values
#         Sharpness=1,
#         kernel='Gaussian',
#         scaleInterp=1.)
#        
#        Tpoints=opyf.Interpolate.npGrid2TargetPoint2D(grid_xp,grid_yp)
#        Interp=opyf.Interpolate.npInterpolateVTK2D(pointsp[:,0:2],vel[:,0:2],Tpoints,ParametreInterpolatorVTK=interp_params)
#        NANmask=np.ones((h,w))
#        NANmask[np.where(mask==0)]=np.nan
#        UxMPM=opyf.Interpolate.npTargetPoints2Grid2D(Interp[:,0],w,h)*NANmask
#        UyMPM=opyf.Interpolate.npTargetPoints2Grid2D(Interp[:,1],w,h)*NANmask
#        grid_val=np.array((UxMPM**2+UyMPM**2)**0.5)                 
#    #                grid_val = d6py.interpolate2D(pointsp,vel,,paramInterpMPM,dConfig,mask=mask)
#    #                write_csvField(grid_val,writingFieldMPMPath+'/grid_val_'+format(ifile,'04.0f')+'.csv')
#                    
#            #Save Grid_val in a CSV file
#    #            write_csvField(grid_val,writingFieldMPMPath+'/grid_val_'+format(frame_idx-framedeb,'04.0f')+'.csv')
#               
##        plt.pause(0.01)
#                  
#        Hfig=8
# 
#        cap.set(cv2.CAP_PROP_POS_FRAMES, vec1[ifile])
#        ret, frame = cap.read()
#    
#        frame=frame[ROI[1]:(ROI[3]+ROI[1]),ROI[0]:(ROI[2]+ROI[0])]
#        vis=opyf.CLAHEbrightness(frame,70)                
#        
#
#        #Show the image with the extent
#    
#        fig,ax=opyf.opyffigureandaxes(extent=extentr,Hfig=4, unit='m',num='a')
#        ax.imshow(cv2.cvtColor(vis, cv2.COLOR_BGR2RGB),extent=extentr)
#                       
#        ax.plot(vecX,H,'--k',linewidth=3,label='H(x) - MPM Sand 6 simulation')
#        ax.plot(coordinatesDP[1][1],variablesDP[0][1][ifile,:],'-k',linewidth=4,label='H(x) - Experiment',alpha=0.5)
#        ax.plot([posDoor[ifile,0],posDoor[ifile,0]],[posDoor[ifile,1]-0.45*dConfig['box'][1]-0.007,posDoor[ifile,1]+0.45*dConfig['box'][1]],'-')
#                    #set the ROI limits      
#        ax.set_xlim([extentr[0],extentr[1]])
#        ax.set_ylim([extentr[2],extentr[3]]) 
##        plt.pause(0.1)
#    
#        plt.savefig(ImgCheckDepth+'/'+format(ifile,'04.0f')+'.png')
#        d6Depth.append(H)
#        ifile+=1
#    writingModDepthPath=out_Opyf+listOutOpyf[j]+'/DATA_sand6'
#    opyf.mkdir2(writingModDepthPath)
#    if door=='with':
#        outputname='/Run_'+format(j,'02.0f')+'_3D_Door_muRigid='+str(muRigid)+'_delta_mu='+format(delta_mu,'.3f')+'TypicalLength='+format(TypicalLength,'.3f')+'_model_depths'+'_substeps_'+str(substeps)
#    else:
#        outputname='/Run_'+format(j,'02.0f')+'_3D_no_Door_start_at_0.13_s_delta_mu='+format(delta_mu,'.3f')+'TypicalLength='+format(TypicalLength,'.3f')+'_model_depths'+'_substeps_'+str(substeps)
#
#    filename=writingModDepthPath+outputname+'.hdf5'
#    d6Depth=np.array(d6Depth)
#    opyf.hdf5_Write(filename,[['T [s]',np.arange(ifile)/dConfig['fps']],['X [m]',vecX]],[['Depth [m]',d6Depth]])
#    copyfile(configFileout, writingModDepthPath+outputname+'.conf')
    
    
    
#%% Write the field of solid fraction and velocity field and stress field




    
#%%                        
     

#            scalebar = ScaleBar(1,length_fraction=0.3)
#            scalebar.box_alpha=0.
#            plt.gca().add_artist(scalebar)
#
#               
#
#                Time=float((frame_idx-frameref))/sdictE['fps']  
#  
#                dech,decv=-0.005,0.005
#                ax.text(0.1+dech,0.08+decv,  'Type : '+str(sdictE['grainType'])+ r'~~Time: ' + format(Time,'1.2f') +' s')
#                l=ax.legend(bbox_to_anchor=((0.09+dech-extentr[0])/(extentr[1]-extentr[0]),(0.068-extentr[2])/(extentr[3]-extentr[2])),loc=6,fontsize=10)
#                l.get_frame().set_linewidth(0.0)
#                l.get_frame().fill=False
##                ax.text(0.1+dech,0.04+decv,)
#                ax.text(0.1+dech,0.04+decv,r'$\mu$ = ' + format(sdictE['mu'],'1.2f') + r'$~~~\Delta \mu$ = ' + format(delta_mu,'1.2f') + r'$~~~Slope$ = '+ format(sdictE['Slope'],'1.0f') +'$^\circ$')
##                ax.text(0.22+dech,0.04+decv,r'$Slope$ ='+ format(sdictE['Slope'],'1.0f') +'$^\circ$')
#                      
##                ax.text(0.23+dech,0.08+decv,r'Time: ' + format(Time,'1.2f') )
#                     
#                     
#                
#        #        fig.set_size_inches(4.5, 3)   
#                
#                plt.savefig(writingImgPath2+'/Step'+format(frame_idx-framedeb,'04.0f')+'.png',format='png',dpi=200) 
#                
##%
#
#                fig2=plt.figure('fig',figsize=(Lfig, Hfig))
#                fig2.clf()
#                ax = plt.Axes(fig2, [0.1, 0.25, 0.85, 0.75])
#                fig2.add_axes(ax)
#                extent=[np.min(grid_xp),np.max(grid_xp),np.min(grid_yp),np.max(grid_yp)]
#                infoPlot={'cmap' : cmap,
#                          'markersize' : 0.3,
#                          'contourlim' : contourlim,
#                          'vlim' : vlim,
#                          'label' : 'Sand6 - MPM - Simulation',
#                          'strpoints': 'Count MPM particles ',
#                          'ifile':frame_idx,
#                          'DP':{'vecX' : vecX, 'H' : H},
#                          'extentr':extentr,
#                          'axes':None}  
#        
#                d6py.d6imshow(grid_val, grid_xp, grid_yp,pointsp,extent,dConfig=dConfig,fig=fig2,infoPlot=infoPlot,sdictE=sdictE)
# 
#                plt.savefig(writingImgPath+'/Step'+format(frame_idx-framedeb,'04.0f')+'.png',format='png',dpi=200) 
#
#                plt.pause(0.01)
#
#    
    
    






