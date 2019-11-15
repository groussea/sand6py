#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 10:48:34 2019

@author: Gauthier
"""
import numpy as np
import csv
#%%
def findHzeroAndMeanVel(Field,Ux,Uy,X,Y,vecXexpD,expD_in):
    import cv2
    h,w=Field.shape
    ret, thresh = cv2.threshold(Field, 0.05, 255, 0)
    cv2.dilate(thresh,kernel=np.ones((3,3)))
    thresh=np.array( thresh,dtype='uint8') 
    contours, hie= cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE )

    drawing = np.zeros((thresh.shape[0], thresh.shape[1]), dtype=np.uint8)
    if len(contours)>0:
        drawing=cv2.drawContours(drawing, contours, np.argmax([len(i) for i in contours]), (255,255,255), 1, cv2.LINE_8, hie, 0)

    Hzero=np.zeros(w)
    indHfromXtovecXexpD=np.zeros(w,'uint16')
    indHfromYtoexpD=np.zeros(w,'uint16')
    deltaH=np.zeros(w)
    Uxmean,Uymean=np.zeros(w),np.zeros(w)
    for j in range (0,len(drawing[0,:])):
        indHfromXtovecXexpD[j]=np.argmin((X[j]-vecXexpD)**2)
    for i in range (0,len(drawing[:,0])):
        indHfromYtoexpD[j]=np.argmin((Y[i]-expD_in)**2)
    for j in range (0,len(drawing[0,:])):
        for i in np.arange(len(drawing[:,0])-1,-1,-1):
            if drawing[i,j] !=0 :
                Hzero[j]=Y[i]
                Uxmean[j]=np.nanmean(Ux[indHfromYtoexpD[i]:i,j])
                Uymean[j]=np.nanmean(Uy[indHfromYtoexpD[i]:i,j])
                break
        if Hzero[j]==0:
            np.argmin(X[j]-vecXexpD)
            Hzero[j]=expD_in[indHfromXtovecXexpD[j]]
        deltaH[j]=expD_in[indHfromXtovecXexpD[j]]-Hzero[j]
        
    return Hzero,deltaH,Uxmean,Uymean

def extractInfoVTK(filename):  
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.Update()
    data=reader.GetOutput()           
    pdata=data.GetPointData()
    #Convert volocities into Numpy Array
    velname=pdata.GetArrayName(1)
    velarray=pdata.GetArray(velname)
    npvelArr=vtk_to_numpy(velarray) 
    #Convert volumes into Numpy Array
    volname=pdata.GetArrayName(0)
    npvolArr=vtk_to_numpy(pdata.GetArray(volname)) 
    #Convert points into Numpy Array
    pointsArr=data.GetPoints().GetData()
    nppointsArr=vtk_to_numpy(pointsArr) 
    
    return npvelArr, npvolArr, nppointsArr

def extractInfoVTKprimalfields(filename):  
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.Update()
    data=reader.GetOutput()           
    pdata=data.GetPointData()
    #Convert volocities into Numpy Array
    velname=pdata.GetArrayName(1)
    velarray=pdata.GetArray(velname)
    npvelArr=vtk_to_numpy(velarray) 
    #Convert solid fraction into Numpy Array
    volname=pdata.GetArrayName(0)
    npPhiArr=vtk_to_numpy(pdata.GetArray(volname)) 
    #Convert stresses into Numpy Array
    strressname=pdata.GetArrayName(2)
    npStressesArr=vtk_to_numpy(pdata.GetArray(strressname))     
    
    #Convert points into Numpy Array
    pointsArr=data.GetPoints().GetData()
    nppointsArr=vtk_to_numpy(pointsArr) 
    
    return npvelArr, npPhiArr, npStressesArr, nppointsArr


def extractInfoVTKprimalfields3D(filename):  
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.Update()
    data=reader.GetOutput()           
    pdata=data.GetPointData()
    #Convert volocities into Numpy Array
    velname=pdata.GetArrayName(1)
    velarray=pdata.GetArray(velname)
    npvelArr=vtk_to_numpy(velarray) 
    #Convert solid fraction into Numpy Array
    phiname=pdata.GetArrayName(0)
    npPhiArr=vtk_to_numpy(pdata.GetArray(phiname)) 
    #Convert stresses into Numpy Array
 
    #Convert points into Numpy Array
    pointsArr=data.GetPoints().GetData()
    nppointsArr=vtk_to_numpy(pointsArr) 
    
    return npvelArr, npPhiArr, nppointsArr


def readConfigFile(configFile):

    f = open(configFile,'r')
    Conf=[]
    for row in f: 
        temp=row.split()
        Conf.append(temp)
    #Extract the different informations 
    f.close()
    Conf=np.array(Conf)
      
#Exemple of Conf file 
# ./apps/d6 1.0-review on release [2017-10-12 14:26]
#fps	150
#substeps	1
#box	0.637	0.135	
#res	60	20	
#nSamples	3
#randomize	0
#volMass	1473
#viscosity	0
#gravity	0	-9.81	
#phiMax	1
#mu	0.43
#delta_mu	0.22
#I0	0.279
#grainDiameter	0.00046
#muRigid	0.5
#cohesion	0
    dictConf={}
    #Create a dictionnary with all the infos we need
    for ic in range(1,len(Conf)):
        if Conf[ic][0]=='box' or Conf[ic][0]=='res' or Conf[ic][0]=='initialOri':
            if len(Conf[ic])==3:
                dictConf[Conf[ic][0]]=[float(Conf[ic][1]),float(Conf[ic][2])]
            else:
                dictConf[Conf[ic][0]]=[float(Conf[ic][1]),float(Conf[ic][2]),float(Conf[ic][3])]
        elif Conf[ic][0]=='scenario' or Conf[ic][0]=='boundary' or Conf[ic][0]=='base_dir': 
            dictConf[Conf[ic][0]]=Conf[ic][1]
        else:    
            dictConf[Conf[ic][0]]=float(Conf[ic][1])
        
            
 #Add to the dictionnay the scale rules for the distances and the velocities                           
    resX=int(dictConf['res'][0])
    resY=int(dictConf['res'][1])
    boxX=float(dictConf['box'][0])
    boxY=float(dictConf['box'][1])
    if len(dictConf['res'])==2:
        gridscale=np.min([boxX/resX,boxY/resY]) #the grid scale is the minimum grid length
 #the velocities are usually scaled with the wave velocity
    else:
        resZ=int(dictConf['res'][2])
        boxZ=float(dictConf['box'][2])
        gridscale=np.min([boxX/resX,boxY/resY,boxZ/resZ])
    
   
    dictConf['gridScale']= gridscale 
    dictConf['velScale']= np.sqrt(gridscale*9.81)     
    return dictConf    

def modifyConfigFile(configFile,newConfigFile,parameter,value):

    f = open(configFile,'r')
    Conf=[]
    for row in f: 
        temp=row.split()
        Conf.append(temp)
    #Extract the different informations 
#    Conf=np.array(Conf)
    f.close()
    
    ofile  = open(newConfigFile, "w")
    writer = csv.writer(ofile, delimiter='\t')
        
    
    for row in Conf:
        a=row[0].find(parameter)
        if (a!=-1) & (len(parameter)==len(row[0])):
            if parameter=='res' or parameter=='nFrames' or parameter=='fps' or parameter=='substeps' or parameter=='nSamples':
                formatstr='.0f'
            else:
                formatstr='.6f'
            print(row)
            if len(value)==1:                  
                row[1]=format(value[0],formatstr)
            elif len(value)==2:
                row[1]=format(value[0],formatstr)
                row[2]=format(value[1],formatstr)
            elif len(value)==3:
                row[1]=format(value[0],formatstr)
                row[2]=format(value[1],formatstr)
                row[3]=format(value[2],formatstr) 
            elif type(value)==str:
                row[1]=value
                del(row[2:])
            print('Changed in '); print(row)
        writer.writerow(row)

    ofile.close()

#    print(configFile+'--- Modified in' +newConfigFile) 
    
    
def RatioAndOrigin(dConfig,dim=3):
    Ratio=(dConfig['columnLength']/(dConfig['box'][0]))
    Ory=float(dConfig['box'][dim-1])*0.1
    Orx=float(dConfig['box'][0])*(Ratio)   
    return Ratio, Orx, Ory
    
def fromVTKtoscaledDataParticles3D(VTKfilename,dConfig):
    Ratio, Orx, Ory=RatioAndOrigin(dConfig)
    
    npvelArr,npvolArr,nppointsArr=extractInfoVTK(VTKfilename)
    vel=npvelArr*dConfig['velScale']
    vol=npvolArr*(dConfig['gridScale'])**3
    points=nppointsArr*dConfig['gridScale']
    
    pointsp=np.zeros_like(points)
    pointsp[:,0]=points[:,0]-Orx
    pointsp[:,2]=points[:,2]-Ory   


    return pointsp,vel,vol        

def fromVTKtoscaledDataParticles2D(VTKfilename,dConfig):
    Ratio, Orx, Ory=RatioAndOrigin(dConfig,dim=2)
    
    npvelArr,npvolArr,nppointsArr=extractInfoVTK(VTKfilename)
    vel=npvelArr*dConfig['velScale']
    vol=npvolArr*(dConfig['gridScale'])**3
    points=nppointsArr*dConfig['gridScale']
    
    pointsp=np.zeros_like(points)
    pointsp[:,0]=points[:,0]-Orx
    pointsp[:,1]=points[:,1]-Ory   


    return pointsp,vel,vol      
    ###translate the origine on the points and the grid (comment if isnot dConfig['columnLength'])

    #% load phi in primal field
def fromVTKtoscaledDataFields3D(VTKfilename,dConfig):        
#    npvelArr, npPhiArr, npStressesArr, nppointsArr=tools.extractInfoVTKprimalfields(VTKfilename)
    Ratio, Orx, Ory=RatioAndOrigin(dConfig)
    npvelArr, npPhiArr, nppointsArr=extractInfoVTKprimalfields3D(VTKfilename)
     
    nppointsArr=nppointsArr*dConfig['gridScale']
    points_grid=np.zeros_like(nppointsArr)             
    points_grid[:,0]=nppointsArr[:,0]-Orx
    points_grid[:,2]=nppointsArr[:,2]-Ory 
    resX=int(dConfig['res'][0])
    resY=int(dConfig['res'][1])
    resZ=int(dConfig['res'][2])
    grid_x=np.reshape(points_grid[:,0],(resX+1,resY+1,resZ+1))
    grid_y=np.reshape(points_grid[:,1],(resX+1,resY+1,resZ+1))
    grid_z=np.reshape(points_grid[:,2],(resX+1,resY+1,resZ+1))
    reshaped_Phi=np.reshape(npPhiArr,(resX+1,resY+1,resZ+1))
    velx=np.reshape(npvelArr[:,0],(resX+1,resY+1,resZ+1))
    vely=np.reshape(npvelArr[:,1],(resX+1,resY+1,resZ+1))
    velz=np.reshape(npvelArr[:,2],(resX+1,resY+1,resZ+1))
    return [grid_x, grid_y, grid_z], [velx,vely,velz],reshaped_Phi


def fromVTKtoscaledDataFields2D(VTKfilename,dConfig):        
#    npvelArr, npPhiArr, npStressesArr, nppointsArr=tools.extractInfoVTKprimalfields(VTKfilename)
    Ratio, Orx, Ory=RatioAndOrigin(dConfig,dim=2)
    npvelArr, npPhiArr, npStressesArr, nppointsArr=extractInfoVTKprimalfields(VTKfilename)
     
    nppointsArr=nppointsArr*dConfig['gridScale']
    points_grid=np.zeros_like(nppointsArr)             
    points_grid[:,0]=nppointsArr[:,0]-Orx
    points_grid[:,1]=nppointsArr[:,1]-Ory 
    resX=int(dConfig['res'][0])
    resY=int(dConfig['res'][1])

    grid_x=np.reshape(points_grid[:,0],(resX+1,resY+1))
    grid_y=np.reshape(points_grid[:,1],(resX+1,resY+1))

    reshaped_Phi=np.reshape(npPhiArr,(resX+1,resY+1))
    velx=np.reshape(npvelArr[:,0],(resX+1,resY+1))
    vely=np.reshape(npvelArr[:,1],(resX+1,resY+1))

    return [grid_x, grid_y], [velx,vely],reshaped_Phi

def csv_write(realpoints,realvel,realvol,csvPath):
    npoints=len(realvel)
    #Write the the scaled values in csv files         
    f=open(csvPath,'w')
    writer = csv.DictWriter(f, fieldnames= ['N','X','Y','Z','Ux','Uy','Uz','Volume'])
    writer.writeheader()        
    for ip in range(0,npoints):
        writer.writerow({'N' : ip,'X' : realpoints[ip][0] ,'Y' : realpoints[ip][1] ,'Z' : realpoints[ip][2] ,'Ux' : realvel[ip][0],'Uy' : realvel[ip][1],'Uz' : realvel[ip][1],'Volume' : realvol[ip]}) 

def csv_read(csvPath):
    reader = csv.reader(open(csvPath,'rb'),csv.QUOTE_NONNUMERIC)
    index=0
    headerDatas=[]
    datas=[]
    for row in reader: 
        if index==0:
            temp=[item for number, item in enumerate(row)]
            headerDatas.append(temp)
        if index>0:
            temp=[float(item) for number, item in enumerate(row)]
            datas.append(temp)
        index+=1
    datas=np.array(datas, dtype=np.float)
    headerDatas=(np.array(headerDatas))
    
    Ux=datas[:,np.where(headerDatas[0] == 'Ux')[0][0]] 
    Uy=datas[:,np.where(headerDatas[0] == 'Uy')[0][0]] 
    Uz=datas[:,np.where(headerDatas[0] == 'Uy')[0][0]] 
    Vel=np.array([Ux,Uy,Uz]).T
    X=datas[:,np.where(headerDatas[0] == 'X')[0][0]]
    Y=datas[:,np.where(headerDatas[0] == 'Y')[0][0]]
    Z=datas[:,np.where(headerDatas[0] == 'Z')[0][0]]
    points=np.array([X,Y,Z]).T
    Volumes=datas[:,np.where(headerDatas[0] == 'Volume')[0][0]]

    return points, Vel, Volumes

#
def readDoorFile(doorF):
    import csv
    headerDatas=[]
    datas=[]
    index=0
    with open(doorF, newline='') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')    
        for row in reader: 
            if index==0:
                temp=[item for number, item in enumerate(row)]
                headerDatas.append(temp)
            if index>0:
                temp=[float(item) for number, item in enumerate(row)]
                datas.append(temp)
            index+=1
    datas=np.array(datas, dtype=np.float)
    headerDatas=(np.array(headerDatas))  
    return headerDatas,datas
  

#%%
#def determineZeroVelocityContour():
    
