# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 16:30:42 2015

@author: thomas
"""
import numpy as np
from scipy import fftpack
from PIL import Image
from PIL import ImageFilter
from PIL import ImageEnhance
from math import pi

def Nearest2(s):
    return 2**((s-1).bit_length()+1)
    
def FigOfMerit(Io,If,sx,sy,Mmask):
    If = If*Mmask
    Io = Io*Mmask
    
    N = np.sum(np.sum(Mmask,1),0)
    If = If/np.sum(np.sum(If,1),0)
    Io = Io/np.sum(np.sum(Io,1),0)
    Ion = np.sum(np.sum(Io,1),0)
    
    eta = 0.0
#    
#    for i in xrange(0,sx):1e-8
#        for j in xrange(0,sy):
#            if(Io[j,i] != 0.0):
#                eta = eta + (If[j,i]/Io[j,i])

    eta = ((If-Io)**2)/(Ion**2)
    eta = np.sum(np.sum(np.sqrt(eta/N),1),0)      
   
    return eta
    
def gaussian_curvature(If,Io):
    Zy, Zx = np.gradient(If-Io)                                                     
    Zxy, Zxx = np.gradient(Zx)                                                  
    Zyy, _ = np.gradient(Zy)                                                    
    K = (Zxx*Zyy-(Zxy** 2))/(1+(Zx**2)+(Zy**2))**2             
    return K
    
def RMS_Error(Io,If,Mmask):
    error = 2
    
    return error


def Foward_Propagation(field_in, lens_phase, xi, d, f, wavelength):
    k = 2*pi/wavelength
    
    phi_1 = np.fft.fftshift(np.fft.fft2(field_in))*np.exp(1j*k*d*np.sqrt(1-wavelength*xi))
    phi_2 = np.convolve(phi_1,lens_phase)
    phi_3 = phi_2*np.exp(1j*k*f*np.sqrt(1-wavelength*xi))
    field_out = np.fft.ifft(np.fft.ifftshift(phi_3))
    return field_out

def Lens_Phase(wavelength, f, r):
    k = 2*pi/wavelength    
    lens = np.exp(1j*k*(r**2)/(2*f))
    lens_phase = np.fft.fftshift(np.fft.fft2(lens))
    return lens_phase
    
def Backwards_Propagation(field_out, lens_phase, xi, d, f, wavelength):
    k = 2*pi/wavelength
    phi_3 = np.fft.shift(np.fft.fft2(field_out))
    phi_2 = phi_3*np.exp(-1j*k*f*np.sqrt(1-wavelength*xi))
    phi_1 = np.convolve(lens_phase,phi_2)
    field_in = np.fft.ifft(np.fft.ifftshift(phi_1))*np.exp(-1j*k*d*np.sqrt(1-wavelength*xi))
    return field_in

m = 0.51
alpha = 0.05
R = 0.0001
L = 160
B = 100
mu = 0

eta = 1000
errors = np.empty(0)
sigX = 80.0
sigY = 80.0
count = 1

name = "Vertical/10_circ"
im = Image.open(name+".bmp").convert('LA')
data = np.array(im.getdata())[:,0].astype(float)
sx,sy = im.size
data = np.reshape(data,(sy,sx))

name = name+"_holo.bmp"

blur = im.filter(ImageFilter.GaussianBlur(50))
blurdata = np.array(blur.getdata())[:,0].astype(float)
blurdata = np.reshape(blurdata,(sy,sx))

#nx = Nearest2(sx)
#ny = Nearest2(sy)
nx = sx    
ny = sy
xGrid = np.arange(1,nx+1,1)
yGrid = np.arange(1,ny+1,1)

Io = np.zeros((ny,nx))
Io[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2] = data

"""OMARF"""
#Iosqrt = np.sqrt(data**2 + 20**2)
#Iosqrt = np.sqrt((Iosqrt/np.max(Iosqrt))*255)
"""IMAGE NOW OFFSET"""
Iosqrt = np.sqrt(Io)

Smask = np.zeros((ny,nx))
Nmask = np.zeros((ny,nx))
Mmask = np.zeros((sx,sy))

temp = (blurdata != 0).astype(int)
Smask[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2] = temp
temp = (Smask == 0).astype(int)
Nmask = temp
temp = (data != 0).astype(int)
Mmask = temp
temp = None

#Smask = np.zeros((ny,nx))
#Smask[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2] = 1
#Nmask = np.ones((ny,nx))>>> >>> sinImg.save('sin_x_100_0.png')
#Nmask[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2] = 0

#Ko = pi*np.random.normal(0,pi/2,(ny,nx))
Ko = pi*np.random.rand(ny,nx)
#Ko = ((Ko)%(2*pi)).astype(complex)

#Ko = 4*R*((alpha*(xGrid-nx/2)**2)[np.newaxis,:]+(1-alpha)*((yGrid-ny/2)**2)[:,np.newaxis])
#Ko = ((Ko)%(2*pi))
#Ko = Ko + L*((xGrid-nx/2)[np.newaxis,:]*np.cos(mu)+(yGrid-ny/2)[:,np.newaxis]*np.sin(mu))
#K2 = B*(((xGrid-nx/2)**(2.0))[np.newaxis,:]+((yGrid-ny/2)**(2.0))[:,np.newaxis])**(0.5)
#Ko = ((Ko)%(2*pi))

Ao = np.exp(-((xGrid-nx/2.0)**2.0)[np.newaxis,:]/(2.0*sigX**2.0)-((yGrid-ny/2.0)**2.0)[:,np.newaxis]/(2.0*sigY**2.0))
#Ao = np.sqrt(Ao.astype(np.uint8))
#Ao = np.zeros((ny,nx))
#for i in xrange(0,np.size(xGrid)):
#    for j in xrange(0, np.size(yGrid)):
#        if (xGrid[i]-nx/2.0 < 3*sigX and yGrid[j]-ny/2.0 < 3*sigY):
#            Ao[j,i] = np.exp(-((xGrid[i]-nx/2.0)**2.0)/(2.0*sigX**2.0)-((yGrid[j]-ny/2.0)**2.0)/(2.0*sigY**2.0))
#        else:
#            Ao[j,i] = 0
Ao = np.sqrt(Ao*255)
#Ao[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2] = Gaussian[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2]

Ein = np.sqrt(np.abs(Ao))*np.exp(1j*Ko)
Eout = np.fft.fftshift(np.fft.fft2(Ein))
G = (m*Iosqrt*Smask+(1-m)*np.abs(Eout)*Nmask)*np.exp(1j*np.angle(Eout))
g = np.fft.ifft2(np.fft.ifftshift(G))    
Ein = np.abs(Ao)*np.exp(1j*np.angle(g))

print("Starting GS Portion")
while (eta>0.00001 and count<50):
    Eout = np.fft.fftshift(np.fft.fft2(Ein))
    If = np.abs(Eout)**2 
    If = (If/np.max(If))*255    
    G = Iosqrt*np.exp(1j*np.angle(Eout))
    g = np.fft.ifft2(np.fft.ifftshift(G))    
    Ein = np.abs(Ao)*np.exp(1j*np.angle(g))
    eta = FigOfMerit(If[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],Io[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],sx,sy,Mmask)     
    print("error = " + str(eta) +"  after "+str(count)+" itterations.")
    count = count +1
    if (count%10 == 0):
        errors = np.append(errors,eta)


print("Starting first Signal Region")
while (eta>0.0000001 and count<00):
    Eout = np.fft.fftshift(np.fft.fft2(Ein))
    If = np.abs(Eout)**2 
    If = (If/np.max(If))*255  
    G = (m*Iosqrt*Smask+(1-m)*np.sqrt(If)*Nmask)*np.exp(1j*np.angle(Eout))
    g = np.fft.ifft2(np.fft.ifftshift(G))    
    Ein = np.abs(Ao)*np.exp(1j*np.angle(g))
    eta = FigOfMerit(If[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],Io[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],sx,sy,Mmask)     
    print("error = " + str(eta) +"  after "+str(count)+" itterations.")
    count = count +1
    if (count%10 == 0):
        errors = np.append(errors,eta)

#Smask = np.zeros((ny,nx))Inbox
#Nmask = np.zeros((ny,nx))
#temp = (blurdata != 0).astype(int)
#Smask[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2] = temp
#temp = (Smask == 0).astype(int)
#Nmask = temp
#temp = None

#KinoR1 = np.angle(Ein)[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2]
#Kino_imgR1 = ((KinoR1-np.min(KinoR1))/(np.max(KinoR1)-np.min(KinoR1)))*255
#phaseImR1 = Image.fromarray(Kino_imgR1.astype(np.uint8))
#RimgR1 = If
#RimgR1Scaled = (((RimgR1-np.min(RimgR1))/(np.max(RimgR1)-np.min(RimgR1)))*255).astype(np.uint8)


#Smask[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2] = 1
#Nmask = np.ones((ny,nx))
#Nmask[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2] = 0

#print "Starting second Signal Region"
#while (eta>0.00001 and count<600):
#    Eout = np.fft.fftshift(np.fft.fft2(Ein))
#    If = np.abs(Eout)**2 
#    If = (If/np.max(If))*255
#    #outputPower = m*np.sum(np.sum(np.abs(Eout)**2*Smask,1),0)+(1-m)*np.sum(np.sum(np.abs(Eout)**2*Nmask,1),0)
#    #print "Output power = "+str(outputPower)+"   Max Eout = " + str(np.max(np.abs(Eout)**2)) 
#    #If = If/np.max(If)    
#    G = (m*Iosqrt*Smask+(1-m)*np.sqrt(If)*Nmask)*np.exp(1j*np.angle(Eout))
#    #G = Iosqrt*np.exp(1j*np.angle(Eout))
#    g = np.fft.ifft2(np.fft.ifftshift(G))
#    Ein = np.abs(Ao)*np.exp(1j*np.angle(g))    
#    #Ein = np.abs(Ao)*np.exp(1j*(np.angle(g)%(1.4*pi)))
#    eta = FigOfMerit(If[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],Io[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],sx,sy,Mmask)     
#    print "error = " + str(eta) +"   after "+str(count)+" itterations."
#    count = count +1

print("MRAF finished after "+str(count)+" iterations.")

RimgR2 = If
RimgR2 = (((RimgR2-np.min(RimgR2))/(np.max(RimgR2)-np.min(RimgR2)))*255).astype(np.uint8)


Rimg = If[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2]
Rimg = (((Rimg-np.min(Rimg))/(np.max(Rimg)-np.min(Rimg)))*255).astype(np.uint8)

roughness = np.sum(np.sum(gaussian_curvature(Rimg,data)**(2.0),1),0)/float((sx*sy)) 
print( "Roughness of intensities = "+str(roughness)+" .")  
Kino = np.angle(Ein)[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2]
#Kino = (Kino)%(1.2*pi)

Kino_img = ((Kino-np.min(Kino))/(np.max(Kino)-np.min(Kino)))*255.0
phaseIm = Image.fromarray(Kino_img.astype(np.uint8))
#contrast_phaseIm = Image.fromarray(Kino_img.astype(np.uint8))
#cont = ImageEnhance.Contrast(contrast_phaseIm)
#cont.enhance(0.3)

phaseIm.save(name+'.bmp')

#Kino_mod = (Kino+pi)*(1.4/2.0)
#Final_Check = np.fft.fftshift(np.fft.fft2(Ao*np.exp(1j*Kino_mod)))
#Final_Check = np.abs(Final_Check)




