#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:04:39 2019

@author: Gauthier Rousseau
"""

import numpy as np
import cv2

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


def mask_collapses(frame,OR,sdictE): 

    h, w = frame.shape[:2] 
   
#      Process
    if len(np.shape(frame))==3:           
        frame_gray=cv2.cvtColor(frame,cv2.COLOR_BGR2GRAY)
    else:
        frame_gray=np.copy(frame)
#    frame_gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    th, BW = cv2.threshold(frame_gray,sdictE['maskTh']*255,255,cv2.THRESH_BINARY)
    BW=cv2.bitwise_not(BW)  
    kernel = np.ones((3,3),np.uint8)
    BW=cv2.morphologyEx(BW, cv2.MORPH_OPEN, kernel)
       #fill the holes
    im_floodfillL=BW.copy()
       # Floodfill from point (0, 0)
    mask = np.zeros((h+2, w+2), np.uint8)
    cv2.floodFill(im_floodfillL, mask, (0,0), 255)

        #Invert floodfilled image
    im_floodfillL=cv2.bitwise_not(im_floodfillL)
    
    im_floodfillR=BW.copy()
       # Floodfill from point (0, 0)
    mask = np.zeros((h+2, w+2), np.uint8)
    cv2.floodFill(im_floodfillR, mask, (w-1,0), 255)

        #Invert floodfilled image
    im_floodfillR=cv2.bitwise_not(im_floodfillR)
            
    BW= BW |  (im_floodfillL & im_floodfillR)
        #fill isolated points in O
    im_floodfill2=cv2.bitwise_not(BW)   
    mask = np.zeros((h+2, w+2), np.uint8)
    cv2.floodFill(im_floodfill2, mask, (w-1,h-1), 255)   
#    im_floodfill_inv2=cv2.bitwise_not(im_floodfill2)
    mask = np.zeros((h+2, w+2), np.uint8)
    BW= BW & im_floodfill2    
    BW = cv2.dilate(BW,kernel,iterations = 2)

#    frame_gray[:,:]=frame_gray[:,:]*(BW[:,:]/255)
#   affichage du profil

    H=np.zeros(w)
    vecX=np.zeros(w)
    for j in range (0,len(BW[0,:])):
        for i in range (0,len(BW[:,0])):
            vecX[j]=(j-OR[0])*sdictE['scale']
            if BW[i,j] !=0 :
                H[j]=(-i+OR[1])*sdictE['scale']
                break

    # coorcetion on the vertical origin from the avrage H on the 
    return BW, vecX, H

