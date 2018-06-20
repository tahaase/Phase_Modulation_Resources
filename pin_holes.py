# -*- coding: utf-8 -*-
"""
Created on Thu May 05 16:48:46 2016

@author: thaa191
"""
import numpy as np
from PIL import Image
from scipy import signal

def Blazed_Grating_H(Period,sX,sY):
    t = np.linspace(0,sX,sX)
    saw = signal.sawtooth(2*np.pi*Period*t)
    grating = np.array([saw,]*sY)
    grating = (grating+1)/2
    return grating
    
def Blazed_Grating_V(Period,sX,sY):
    t = np.linspace(0,sY,sY)
    saw = signal.sawtooth(2*np.pi*Period*t)
    grating = np.transpose(np.array([saw,]*sX))
    grating = (grating+1)/2
    return grating

sX = 792
sY = 600
offset_y = 80
offset_x = 0
#data = np.fromfile('prism_60_dat.dat', dtype = int, sep='\t')
#data = np.reshape(data,(sY,sX))
intDis = np.zeros((sY,sX))


size = 150

for i in range(sX):
    for j in range(sY):
        if(((i-sX/2)**2+(j-sY/2)**2) <= (size**2)):
            intDis[j,i] = 1
            
intDis = intDis
grating_h = Blazed_Grating_H(0,sX,sY)
grating_v = Blazed_Grating_V(offset_y,sX,sY)
total_grating = (grating_h+grating_v)%1
intDis = intDis*total_grating*255
img = Image.fromarray(intDis.astype(np.uint8))

img.save(str(offset_y)+"_pin.bmp")
