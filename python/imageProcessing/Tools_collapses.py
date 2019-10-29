#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:04:39 2019

@author: Gauthier Rousseau
"""

import numpy as np
import cv2

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

