#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 10:48:34 2019

@author: Gauthier Rousseau
"""
import os
import numpy as np
import json
import csv
import vtk
from vtk.util.numpy_support import vtk_to_numpy



#%%

def mkdir2(path):
    if not os.path.isdir(path):
        os.mkdir(path)

def extractInfoVTK(filename):  
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
    reader = vtk.vtkDataSetReader()
    #required to read all file datas
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllTensorsOn()
    reader.SetFileName(filename)
    reader.Update()
    data=reader.GetOutput()           
    pdata=data.GetPointData()
    N_array=pdata.GetNumberOfArrays()
    #Convert points into Numpy Array
    np_points=vtk_to_numpy(data.GetPoints().GetData()) 
    
    #Convert solid fraction into Numpy Array
    np_phi=vtk_to_numpy(pdata.GetArray(pdata.GetArrayName(0))) 
    #Convert volocities into Numpy Array
    np_vel=vtk_to_numpy(pdata.GetArray(pdata.GetArrayName(1))) 

        #Convert d_phi into Numpy Array
    np_d_phi=vtk_to_numpy(pdata.GetArray(pdata.GetArrayName(2)))   

    if N_array>3:
                #convert forces
        np_forces=vtk_to_numpy(pdata.GetArray(pdata.GetArrayName(3))) 
        #Convert stresses into Numpy Array
        np_stresses=vtk_to_numpy(pdata.GetArray(pdata.GetArrayName(4)))     
    
    else:
        # vtk file is in mode 'extractAllField=False'
        np_stresses=None        
        np_forces=None
    
    return  np_points, np_phi, np_vel,  np_stresses, np_d_phi, np_forces, N_array




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
        if Conf[ic][0]=='box' or Conf[ic][0]=='res' or Conf[ic][0]=='initialOri' or Conf[ic][0]=='gravity':
            if len(Conf[ic])==3:
                dictConf[Conf[ic][0]]=[float(Conf[ic][1]),float(Conf[ic][2])]
            else:
                dictConf[Conf[ic][0]]=[float(Conf[ic][1]),float(Conf[ic][2]),float(Conf[ic][3])]
        elif Conf[ic][0]=='scenario' or Conf[ic][0]=='boundary' or Conf[ic][0]=='base_dir': 
            dictConf[Conf[ic][0]]=Conf[ic][1]
            if Conf[ic][0]=='scenario':
                #load scenario options
           
                for op in Conf[ic][2:]:
                    split_op=op.split(':')
                    try:
                        dictConf[split_op[0]]=float(split_op[1])
                    except: 
                        dictConf[split_op[0]]=split_op[1]
        else:    
            dictConf[Conf[ic][0]]=float(Conf[ic][1])

            
 #Add to the dictionnay the scale rules for the distances and the velocities                           
    resX=int(dictConf['res'][0])
    resY=int(dictConf['res'][1])
    boxX=float(dictConf['box'][0])
    boxY=float(dictConf['box'][1])
    if len(dictConf['res'])==2:
        typ_L=np.min([boxX/resX,boxY/resY]) #the grid scale is the minimum grid length
 #the velocities are usually scaled with the wave velocity
    else:
        resZ=int(dictConf['res'][2])
        boxZ=float(dictConf['box'][2])
        typ_L=np.min([boxX/resX,boxY/resY,boxZ/resZ])
    
    typ_G=np.linalg.norm(dictConf['gravity'])
    typ_R=dictConf['volMass']

    dictConf['gridScale']= typ_L
    dictConf['velScale']= np.sqrt(typ_L*typ_G) 
    dictConf['pressureScale']= typ_G*typ_R*typ_L
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
    try:
        Ratio=(dConfig['columnLength']/(dConfig['box'][0]))
    except:
        Ratio=(dConfig['column_length']/(dConfig['box'][0]))
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

def fromVTKtoscaledDataFields3D(VTKfilename,dConfig):        

    Ratio, Orx, Ory=RatioAndOrigin(dConfig)
    np_points, np_phi, np_vel,  np_stresses, np_d_phi, np_forces, N_array =extractInfoVTKprimalfields(VTKfilename)
    
    np_points=np_points*dConfig['gridScale']
    points_grid=np.zeros_like(np_points)             
    points_grid[:,0]=np_points[:,0]-Orx
    points_grid[:,1]=np_points[:,1]-Ory 
    resX=int(dConfig['res'][0])
    resY=int(dConfig['res'][1])
    resZ=int(dConfig['res'][2])
    grid_x=np.reshape(points_grid[:,0],(resX+1,resY+1,resZ+1))
    grid_y=np.reshape(points_grid[:,1],(resX+1,resY+1,resZ+1))
    grid_z=np.reshape(points_grid[:,2],(resX+1,resY+1,resZ+1))

    reshaped_Phi=np.reshape(np_phi,(resX+1,resY+1,resZ+1))
    velx=np.reshape(np_vel[:,0],(resX+1,resY+1,resZ+1))*dConfig['velScale']
    vely=np.reshape(np_vel[:,1],(resX+1,resY+1,resZ+1))*dConfig['velScale']
    velz=np.reshape(np_vel[:,2],(resX+1,resY+1,resZ+1))*dConfig['velScale']
    np_d_phi_reshaped=np.reshape(np_d_phi,(resX+1,resY+1,resZ+1,3))
    
    if N_array>3:
        np_stresses_reshaped=np.reshape(np_stresses,(resX+1,resY+1,resZ+1,9))*dConfig['pressureScale']
        np_forces_reshaped=np.reshape(np_forces,(resX+1,resY+1,resZ+1,3))*dConfig['pressureScale']
    else:
        np_stresses_reshaped=None
        np_forces_reshaped=None

    return [grid_x, grid_y, grid_z], [velx,vely,velz], reshaped_Phi, np_d_phi_reshaped, np_stresses_reshaped, np_forces_reshaped



def fromVTKtoscaledDataFields2D(VTKfilename,dConfig):        
#    npvelArr, npPhiArr, npStressesArr, nppointsArr=tools.extractInfoVTKprimalfields(VTKfilename)
    Ratio, Orx, Ory=RatioAndOrigin(dConfig,dim=2)
    np_points, np_phi, np_vel,  np_stresses, np_d_phi, np_forces, N_array =extractInfoVTKprimalfields(VTKfilename)
    
    np_points=np_points*dConfig['gridScale']
    points_grid=np.zeros_like(np_points)             
    points_grid[:,0]=np_points[:,0]-Orx
    points_grid[:,1]=np_points[:,1]-Ory 
    resX=int(dConfig['res'][0])
    resY=int(dConfig['res'][1])

    grid_x=np.reshape(points_grid[:,0],(resX+1,resY+1))
    grid_y=np.reshape(points_grid[:,1],(resX+1,resY+1))
    

    reshaped_Phi=np.reshape(np_phi,(resX+1,resY+1))
    velx=np.reshape(np_vel[:,0],(resX+1,resY+1))*dConfig['velScale']
    vely=np.reshape(np_vel[:,1],(resX+1,resY+1))*dConfig['velScale']

    np_d_phi_reshaped=np.reshape(np_d_phi,(resX+1,resY+1,3))

    if N_array>3:
        np_stresses_reshaped=np.reshape(np_stresses,(resX+1,resY+1,9))*dConfig['pressureScale']
        np_forces_reshaped=np.reshape(np_forces,(resX+1,resY+1,3))*dConfig['pressureScale']
    else:
        np_stresses_reshaped=None
        np_forces_reshaped=None

    return [grid_x, grid_y], [velx,vely], reshaped_Phi, np_d_phi_reshaped, np_stresses_reshaped, np_forces_reshaped

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

def findFoldersInFolder(main):   
   listFolder,paths=[],[]
   for elem in os.listdir(main):
       if os.path.isdir(main+'/'+elem):
           paths.append(main+'/'+elem)
           listFolder.append(elem)
   return listFolder, paths


def findSubFolders(main,N):
    childs=[]
    listFolder,pathsAll=findFoldersInFolder(main)
    for l,p in zip(listFolder,pathsAll):
        childs.append([l])     
    
    n=0
    subFoldersList=childs.copy()
    
    while n<N:
        childsTemp=[]
        for c in childs:
            listFolder,paths=findFoldersInFolder(main+''.join('/'+str(e) for e in c))
            for l,p in zip(listFolder,paths):
                childsTemp.append(c+[l])
                subFoldersList.append(c+[l])
                pathsAll.append(p)
        n+=1
        childs=childsTemp.copy()
    return subFoldersList, pathsAll

def whereFolder(name,N,main='./'):
    mark=[]
    subFoldersList, pathsAll=findSubFolders(main,N)
    for sub in subFoldersList:
        if sub[-1]==name: mark.append(1)
        else: mark.append(0)
    return mark, subFoldersList, pathsAll


def findOutSand6Paths(main,N):
    childs=[['']]

    n=0
    subFoldersList=childs.copy()
    outPaths,outFolders,pathsAll=[],[],[]
    while n<N:
        childsTemp=[]
        for c in childs:
            listFolder,paths=findFoldersInFolder(main+''.join('/'+str(e) for e in c))
            for l,p in zip(listFolder,paths):

                if os.path.exists(p+'/config') and os.path.exists(p+'/stats.txt'): 
                    outPaths.append(p)
                    outFolders.append(l)
                else:
                    childsTemp.append(c+[l])
                    subFoldersList.append(c+[l])
                    pathsAll.append(p)
                    
        n+=1
        childs=childsTemp.copy()
    listDictConf=[]
    listNumericalRuns=[]
    for f,p in zip(outFolders,outPaths):
        currentNumRun=NumericalRun(p)
        currentNumRun.dConfig['folder']=f
        currentNumRun.dConfig['path']=p
        currentNumRun.dConfig['runNumber']=int(currentNumRun.dConfig['folder'][5])
        listDictConf.append(currentNumRun.dConfig)  
        listNumericalRuns.append(currentNumRun)
    return outPaths,outFolders,listDictConf,listNumericalRuns

def whereSand6OutFromDict(listDictConf,dictConf):
    ind=1
    for d in listDictConf:
        if d==dictConf:
            break
        ind+=1
    return d['path'] 

def whereSand6OutFromFolder(listDictConf,folder):
    ind=1
    for d in listDictConf:
        if folder==d['folder']:
            break
        ind+=1
    return d   

def whereSand6OutFromParms(listRuns,**params):
    listRunsOut=listRuns.copy()
    for p in params:
        selectedRuns=[]
        for d in listRunsOut:
            try:
                if d.dConfig[p]==params[p]:
                    selectedRuns.append(d)
            except:
                print('The parameter is not in config dictionnary')
        listRunsOut=selectedRuns.copy()
    listDictOut=[]
    print('Selected Runs are:')
    for r in listRunsOut:
        print('')
        print('Run ' + r.dConfig['folder'])
        listDictOut.append(r.dConfig)
#        for keys,values in r.dConfig.items():
#            print(keys,':',values)
        
    return listRunsOut, listDictOut


#%%            



class NumericalRun():
    def __init__(self,d6OutFolder,doorFolder=None):

        self.dConfig=readConfigFile(d6OutFolder+'/config')     
        dimSim=len(self.dConfig['res'])
        self.dConfig['dimSim']=len(self.dConfig['res'])
        
        if dimSim==3:
            self.Ratio, self.Orx, self.Ory=RatioAndOrigin(self.dConfig)
            self.Ldoor=float(self.dConfig['box'][2])*0.9
        else:
            self.Ratio, self.Orx, self.Ory=RatioAndOrigin(self.dConfig,dim=2)
            self.Ldoor=float(self.dConfig['box'][1])*0.9
#        if doorFolder is not None:
#            self.headerDatas,self.datas=readDoorFile(doorFolder+'/Run_O'+str(runNumber)+'_door.txt') 
#            self.datas[:,3]=self.datas[:,3]-self.Orx
#            self.datas[:,5]=self.datas[:,5]-self.Ory
        self.d6OutFolder=d6OutFolder
        self.dimSim=dimSim
        self.H=self.Ldoor
        self.vecxMax=np.array([2,2.5,3,3,2.5,3,4,6])*self.H

#        self.runNumber=runNumber
#        self.dConfig['runNunmber']=runNumber
#        self.xmax=self.vecxMax[self.runNumber]
        self.scaleLength=1
        self.strainRateCalculated=False
        self.pressureCalculated=False
        self.velocityCalculated=False
        
    def loadVTK(self, ifile):
        self.ifile=ifile
        self.strainRateCalculated=False
        self.pressureCalculated=False
        self.velocityCalculated=False
        VTKfilename_part = self.d6OutFolder+'/vtk/particles-'+ str(ifile) +'.vtk'  
        VTKfilename= self.d6OutFolder+'/vtk/primal-fields-'+ str(ifile) +'.vtk' 
        if self.dimSim==3: 
           [self.grid_x, self.grid_y, self.grid_z], [self.velx,self.vely,self.velz],self.reshaped_Phi, self.np_d_phi_reshaped, self.np_stresses_reshaped, self.np_forces_reshaped=fromVTKtoscaledDataFields3D(VTKfilename,self.dConfig)   
           self.pointsp,self.vel,self.vol =fromVTKtoscaledDataParticles3D(VTKfilename_part,self.dConfig)

        if self.dimSim==2:
            [self.grid_x, self.grid_y], [self.velx,self.vely],self.reshaped_Phi,  self.np_d_phi_reshaped, self.np_stresses_reshaped, self.np_forces_reshaped=fromVTKtoscaledDataFields2D(VTKfilename,self.dConfig)
            self.pointsp,self.vel,self.vol =fromVTKtoscaledDataParticles2D(VTKfilename_part,self.dConfig)
            self.extentR=np.array([self.grid_x[0,0]-self.dx/2,self.grid_x[-1,0]+self.dx/2,self.grid_y[0,0]-self.dy/2,self.grid_y[0,-1]+self.dy/2])/self.scaleLength
            self.dx=np.absolute(self.grid_x[0,0]-self.grid_x[1,0])
            self.dy=np.absolute(self.grid_y[0,0]-self.grid_y[0,1])
        self.phiT=np.flipud(self.reshaped_Phi.T)    
        



    def scLength(self,length):       
        self.scaleLength=length

    def plotContour(self,ax,**args):
        if self.dimSim==3: 
            self.CS=ax.contour(self.grid_x[:,6,:]/self.scaleLength,self.grid_z[:,6,:]/self.scaleLength,self.reshaped_Phi[:,6,:],**args)
        if self.dimSim==2:
            self.CS=ax.contour(self.grid_x/self.scaleLength,self.grid_y/self.scaleLength,self.reshaped_Phi,**args)
    def calculatePressure(self):
        
        if self.dimSim==2:
            self.P=np.flipud(1/2*(self.np_stresses_reshaped[:,:,0]+self.np_stresses_reshaped[:,:,4]).T)
        if self.dimSim==3:
            self.P=np.flipud(1/3*(self.np_stresses_reshaped[:,:,:,0]+self.np_stresses_reshaped[:,:,:,4]+self.np_stresses_reshaped[:,:,:,8]).T)
        self.pressureCalculated=True
        
    def calculateStrainRate(self):
        
        if self.dimSim==2:
            self.P=np.flipud(1/2*(self.np_stresses_reshaped[:,:,0]+self.np_stresses_reshaped[:,:,4]).T)
        if self.dimSim==3:
            [vxdx,vxdy,vxdz]=np.gradient(self.velx,self.dx,self.dy,self.dz)
            [vydx,vydy,vydz]=np.gradient(self.vely,self.dx,self.dy,self.dz)
            [vzdx,vzdy,vzdz]=np.gradient(self.velz,self.dx,self.dy,self.dz)
            
            (rx,ry,rz)=self.velx.shape
            self.gammaN=np.zeros_like(self.velz)
        
            self.sym=np.array(self.gammaN,dtype=object)
            self.J=np.array(self.gammaN,dtype=object)
            for i in range(rx):
                for j in range(ry):
                   for k in range(rz): 
                    self.J[i,j,k]=np.array([[vxdx[i,j,k], vxdy[i,j,k], vxdz[i,j,k]],[vydx[i,j,k], vydy[i,j,k], vydz[i,j,k]],[vzdx[i,j,k], vzdy[i,j,k], vzdz[i,j,k]]])
                    self.sym[i,j,k]=0.5*(self.J[i,j,k]+self.J[i,j,k].T)
                    self.gammaN[i,j,k]=(0.5**0.5)*np.linalg.norm(self.sym[i,j,k])
                    self.normV=np.linalg.norm([self.velx,self.vely])
        self.strainRateCalculated=True
        
    def calculateNormVelocity(self):
        if self.dimSim==2:
            ax.imshow(self.P,extent=self.extentR,**args)
        if self.dimSim==3:
            self.normV=((self.velx)**2+(self.vely)**2+(self.velz)**2)**0.5
        self.velocityCalculated=True
        
        
    def plotVelocity(self,ax,**args):       
        if self.velocityCalculated==False:
            self.calculateNormVelocity()
        self.IMvel=np.flipud(self.normV.T)
        if self.dimSim==2:
            ax.imshow(self.P,extent=self.extentR,**args)

    def plotStrainRate(self,ax,**args):
        if self.strainRateCalculated==False:
            self.calculateStrainRate()
        self.IMgamma=np.flipud(self.gammaN.T)
        if self.dimSim==2:
            ax.imshow(self.IMgamma*(self.phiT>0.5),extent=self.extentR,**args)
        if self.dimSim==3:
            ax.imshow(self.IMgamma[:,6,:]*(self.phiT[:,6,:]>0.5),extent=self.extentR,**args)        
    
    def plotPressure(self,ax,**args):
        if self.pressureCalculated==False:
            if self.np_stresses_reshaped is not None:
                self.calculatePressure() 
            else:
                print('Warning : it is impossible to calculate the pressure from the sand6 outputs (set exportAllFields in the config file to 1 before the simulation if you wish to export stresses and deduce pressure)')
        if hasattr(self, 'P') :
            if self.dimSim==2:
                ax.imshow(self.P,extent=self.extentR,**args)
            if self.dimSim==3:
                ax.imshow(self.P[:,6,:],extent=self.extentR,**args)
            
    def plotI(self,ax,**args):
        if self.pressureCalculated==False:
            if self.np_stresses_reshaped is not None:
                self.calculatePressure() 
            else:
                print('Warning : it is impossible to calculate the pressure from the sand6 outputs (set exportAllFields in the config file to 1 before the simulation if you wish to export stresses and deduce pressure)')
                
        if self.strainRateCalculated==False:
            self.calculateStrainRate()
            
        if hasattr(self, 'P') and hasattr(self, 'gammaN'):
            if self.dimSim==2:
                self.IMgamma=np.flipud(self.gammaN[:,:].T)
                ax.imshow(self.IMgamma*self.dConfig['grainDiameter']/(self.P/self.dConfig['volMass'])*(self.phiT>0.5),extent=self.extentR,**args)
            if self.dimSim==3:
                self.IMgamma=np.flipud(self.gammaN[:,6,:].T)
                ax.imshow(self.IMgamma*self.dConfig['grainDiameter']/(self.P[:,6,:]/self.dConfig['volMass'])*(self.phiT[:,6,:]>0.5),extent=self.extentR,**args)
            
    
    def plotPoints(self,ax,**args):
        if self.dimSim==3: 
            ax.plot(self.pointsp[:,0]/self.scaleLength,self.pointsp[:,2]/self.scaleLength,**args)
        if self.dimSim==2:
            ax.plot(self.pointsp[:,0]/self.scaleLength,self.pointsp[:,1]/self.scaleLength,**args)



    def plotDoor(self,ax):
        if self.ifile<100:
            X=np.array([self.datas[self.ifile,3],self.datas[self.ifile,3]])/self.scaleLength
            Y=np.array([self.datas[self.ifile,5]+self.Ldoor/2,self.datas[self.ifile,5]-self.Ldoor/2])/self.scaleLength
            ax.plot(X,Y,'k')


    def opyfPointCloudColoredScatter(self,ax,nvec=3000,**args):
       
        from matplotlib.colors import Normalize
        if len(self.pointsp) < nvec:
            N = len(self.pointsp)
        else:
            N = nvec
            print('only '+str(N)+'vectors plotted because length(X) >' + str(nvec))

        ind = np.random.choice(np.arange(len(self.pointsp)), N, replace=False)
        Xc = self.pointsp[ind, :]
        Vc = self.vel[ind, :]
        if self.dimSim==3:
            self.norm_vel=(Vc[:,0]**2+Vc[:,2]**2)**0.5

            norm = Normalize()
            norm.autoscale(self.norm_vel)
            self.im = ax.scatter(Xc[:, 0]/self.scaleLength, Xc[:, 2]/self.scaleLength, c=self.norm_vel,**args)


        
        
      
class ExperimentalRun():

     
     def __init__(self, runNumber=0,mainFolder='.',loadSurfaceElevation=True,loadPoints=False,loadField=False):
     # OPYF LIBRARY IS REQUIRED
        import opyf  
        self.runNumber=runNumber
        
#load json experiment dictionnary        
        JSONpath=mainFolder+'/Granular_Collapses_Experimental_Informations.json'
        in_file = open(JSONpath,"r")
        self.dictExp = json.load(in_file) 
        in_file.close()  
        listExp=np.sort([d for d in self.dictExp ])
        listExpFolders=np.sort([mainFolder+'/'+d for d in self.dictExp ])

        self.dictE=self.dictExp[listExp[runNumber]]
        self.typicalTime=(self.dictE['H']/9.81)**0.5
                
        self.vecxMax=[2,2.5,3,3,2.5,3,4,6,8]
        self.xmax=self.vecxMax[runNumber]

        self.nFrames=int(self.dictExp[listExp[runNumber]]['nFrames'])
        if loadSurfaceElevation==True:
            [[_,self.Time],[_,self.vecXexpD]],[[_,self.expD]]=opyf.hdf5_Read(listExpFolders[runNumber]+'/Run_'+format(runNumber,'02.0f')+'_depths.hdf5')
            self.Time=self.Time-self.Time[0]
        if loadField==True:
            coordinates, variables=opyf.hdf5_Read(listExpFolders[runNumber]+'/Run_'+format(runNumber,'02.0f')+'_rectilinear_velocity_fields.hdf5')        
            self.Time=coordinates[0][1]
            self.X=coordinates[1][1]
            self.Y=coordinates[2][1]
            self.dx=self.dy=self.X[1]-self.X[0]
            self.Ux=variables[0][1]
            self.Uy=variables[1][1]
            self.Ux[np.where(self.Ux==0)]=np.nan
            self.Uy[np.where(self.Uy==0)]=np.nan
            # calculate strain rate
            self.epsilon21=[]
            for k in range(0,len(self.Time)):
                dUxdy=(self.Ux[k][:-1,:]-self.Ux[k][1:,:])/self.dy
                self.dUxdyReshaped=(dUxdy[1:,1:-1]+dUxdy[:-1,1:-1])/2          
                dUydx=(self.Uy[k][:,1:]-self.Uy[k][:,:-1])/self.dx
                self.dUydxReshaped=(dUydx[1:-1,1:]+dUydx[1:-1,:-1])/2    
                self.epsilon21.append(0.5*(self.dUxdyReshaped+self.dUydxReshaped))
            self.epsilon21=np.array(self.epsilon21)

            self.h,self.w=self.Ux[0].shape
            
        self.scaleLength=1.
        

        
     def scLength(self,length=1.):
         #TODO all scalings
         self.scaleLength=length
         
     def plotDepthProfile(self,ax,ifile,**args):
        self.l=ax.plot(self.vecXexpD/self.scaleLength,self.expD[ifile]/self.scaleLength,**args)

     def plotField(self,ax,ifile,Type='velocity_norm',**args):
        Xplot=self.X/self.scaleLength
        Yplot=self.Y/self.scaleLength
        dxp=self.dx/self.scaleLength
        dyp=self.dy/self.scaleLength
        if Type=='velocity_norm':
            self.im=ax.imshow((self.Ux[ifile]**2+self.Uy[ifile]**2)**0.5,extent=[Xplot[0]-dxp/2,Xplot[-1]+dxp/2,Yplot[-1]-dyp/2,Yplot[0]+dyp/2],**args)
        elif Type=='shear_rate':            
            self.im=ax.imshow(self.epsilon21,extent=[Xplot[0]+dxp/2,Xplot[-1]-dxp/2,Yplot[-1]+dyp/2,Yplot[0]-dyp/2],alpha=0.5,**args)
#        elif Type=='inertia':
            

def setFigure(fig,ax,sL,unit='-'):
    ax.set_aspect('auto')
    ax.set_position([0.15,0.15,0.8,0.8]) 
    ax.set_xlabel(r'$x$ ['+ unit + ']' )
    ax.xaxis.set_label_coords(0.5, -0.1)
    ax.set_ylabel(r'$h$  ['+unit + ']' )
    ax.yaxis.set_label_coords(-0.1,0.5 )
    #opyfDisp.ax.xaxis.set_ticklabels([])
#    ax.yaxis.set_ticks([0,0.2,0.4,0.6,0.8,1])
    fig.set_size_inches((8*2,4.5*2))    
    ax.plot([-3/sL,4/sL],[0/sL,0/sL],'k',linewidth=0.5,zorder=1.)
    ax.plot([0,0],[-3/sL,3/sL],'k',linewidth=0.5,zorder=1.) 
    ax.grid(zorder=-0.)
    fig.show()                               

    

#%%
#def determineZeroVelocityContour():
    
