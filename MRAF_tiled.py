# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 12:39:14 2016

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
 
def Zernike_Combination((x,y), m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15):
    
    rho = np.zeros((np.size(y),np.size(x)))    
    phi = np.zeros((np.size(y),np.size(x)))    
    for i in xrange(0,np.size(x)):
        for j in xrange(0,np.size(y)):
            rho[j,i] = (x[i]**2.0+y[j]**2.0)**0.5
            phi[j,i] = np.arctan2(y[j],x[i])
    #Piston
    p1 = 1.0 
    #Tip (x-tilt)
    p2 = rho*np.cos(phi)#2.0*rho*np.cos(phi)
    #Tilt (y-tilt)
    p3 = rho*np.sin(phi)#2.0*rho*np.sin(phi)
    #Defocus
    p4 = (2.0*rho**2.0 - 1.0)#(3.0**0.5)*(2.0*rho**2.0 - 1.0)
    
    #Oblique Astigmatism
    p5 = (rho**2.0)*np.sin(2.0*phi)#(6.0**0.5)*(rho**2.0)*np.sin(2.0*phi)
    
    #Vertical Astigmatism    
    p6 = (rho**2.0)*np.cos(2.0*phi)#(6.0**0.5)*(rho**2.0)*np.cos(2.0*phi)
    
    #Vertical Coma    
    p7 = (3.0*rho**3.0 - 2.0*rho)*np.sin(phi)#(8.0**0.5)*(3.0*rho**3.0 - 2.0*rho)*np.sin(phi)
    
    #Horizontal Coma    
    p8 = (3.0*rho**3.0 - 2.0*rho)*np.cos(phi)#(8.0**0.5)*(3.0*rho**3.0 - 2.0*rho)*np.cos(phi)
    
    #Vertical Trefoil    
    p9 = (rho**3.0)*np.sin(3.0*phi)#(8.0**0.5)*(rho**3.0)*np.sin(3.0*phi)
    
    #Oblique Trefoil
    p10 = (rho**3.0)*np.cos(3.0*phi)#(8.0**0.5)*(rho**3.0)*np.cos(3.0*phi) 
    
    #Primary Spherical    
    p11 = (6.0*rho**4.0 - 6.0*rho**2.0 + 1.0)#(5.0**0.5)*(6.0*rho**4.0 - 6.0*rho**2.0 + 1.0)
    
    #Vertical Secondary Astigmatism    
    p12 = (4*rho**4 - 3*rho**2)*np.cos(2*phi)#(10.0**0.5)*(4*rho**4 - 3*rho**2)*np.cos(2*phi)
    
    #Horizontal Secondary Asstigmatism
    p13 = (4.0*rho**4.0 - 3.0*rho**2.0)*np.sin(2.0*phi)#(10.0**0.5)*(4.0*rho**4.0 - 3.0*rho**2.0)*np.sin(2.0*phi)
    
    #Vertical Quadrafoil    
    p14 = (rho**4.0)*np.cos(4.0*phi)#(10.0*0.5)*(rho**4.0)*np.cos(4.0*phi)
    
    #Oblique Quadrafoil
    p15 = (rho**4.0)*np.sin(4.0*phi)#(10.0*0.5)*(rho**4.0)*np.sin(4.0*phi)
    
    full = m1*p1 + m2*p2 + m3*p3 + m4*p4 + m5*p5 + m6*p6 + m7*p7 + m8*p8 + m9*p9 + m10*p10 + m11*p11 + m12*p12 + m13*p13 + m14*p14 + m15*p15    
    
    return full

   
m = 0.49
mod_val = 2.0
max_val = 255*(mod_val/2.0)
alpha = 0.1
R = 0.0003
L = 136
B = 100

eta = 1000
sigX = 40.0
sigY = 40.0
count = 1

im = Image.open("square_200.png").convert('LA')
data = np.array(im.getdata())[:,0].astype(float)
sx,sy = im.size
data = np.reshape(data,(sy,sx))

blur = im.filter(ImageFilter.GaussianBlur(20))
blurdata = np.array(blur.getdata())[:,0].astype(float)
blurdata = np.reshape(blurdata,(sy,sx))

nx = sx
ny = sy
xGrid = np.arange(1,nx+1,1)
yGrid = np.arange(1,ny+1,1)

Io = np.zeros((ny,nx))
Io[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2] = data
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

Ko = 4*R*((alpha*(xGrid-nx/2)**2)[np.newaxis,:]+(1-alpha)*((yGrid-ny/2)**2)[:,np.newaxis])
Ko = ((Ko)%(2*pi))

Ao = np.exp(-((xGrid-nx/2.0)**2.0)[np.newaxis,:]/(2.0*sigX**2.0)-((yGrid-ny/2.0)**2.0)[:,np.newaxis]/(2.0*sigY**2.0))
Ao = np.sqrt(Ao*255)

Ein = np.sqrt(np.abs(Ao))*np.exp(1j*Ko)
Eout = np.fft.fftshift(np.fft.fft2(Ein))
G = (m*Iosqrt*Smask+(1-m)*np.abs(Eout)*Nmask)*np.exp(1j*np.angle(Eout))
g = np.fft.ifft2(np.fft.ifftshift(G))    
Ein = np.abs(Ao)*np.exp(1j*np.angle(g))

print "Starting GS Portion"
while (eta>0.00001 and count<0):
    Eout = np.fft.fftshift(np.fft.fft2(Ein))
    If = np.abs(Eout)**2 
    If = (If/np.max(If))*255    
    G = Iosqrt*np.exp(1j*np.angle(Eout))
    g = np.fft.ifft2(np.fft.ifftshift(G))    
    Ein = np.abs(Ao)*np.exp(1j*np.angle(g))
    eta = FigOfMerit(If[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],Io[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],sx,sy,Mmask)     
    print "error = " + str(eta) +"  after "+str(count)+" itterations."
    count = count +1

print "Starting first Signal Region"
while (eta>0.0000001 and count<100):
    Eout = np.fft.fftshift(np.fft.fft2(Ein))
    If = np.abs(Eout)**2 
    If = (If/np.max(If))*255  
    G = (m*Iosqrt*Smask+(1-m)*np.sqrt(If)*Nmask)*np.exp(1j*np.angle(Eout))
    g = np.fft.ifft2(np.fft.ifftshift(G))    
    Ein = np.abs(Ao)*np.exp(1j*np.angle(g))
    eta = FigOfMerit(If[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],Io[(ny-sy)/2:(ny+sy)/2,(nx-sx)/2:(nx+sx)/2],sx,sy,Mmask)     
    print "error = " + str(eta) +"  after "+str(count)+" itterations."
    count = count +1

print "MRAF finished after "+str(count)+" iterations."

Rimg = If
Rimg = (((Rimg-np.min(Rimg))/(np.max(Rimg)-np.min(Rimg)))*255).astype(np.uint8)

roughness = np.sum(np.sum(gaussian_curvature(Rimg,data)**(2.0),1),0)/float((sx*sy)) 
print "Roughness of intensities = "+str(roughness)+" ."  

Kino = np.angle(Ein)


""" Add SH correction """

x = np.linspace(-(sx)+1,sx,sx)/sx
y = np.linspace(-(sy)+1,sy,sy)/sy
coeffs = np.array([ 0.054118  ,  0.00397079, -0.00506472, -0.01631324, -0.00201275,
       -0.05382548,  0.00495238,  0.00023831, -0.00322624,  0.00090626,
       -0.00998653,  0.03292146,  0.00012583, -0.01058785, -0.00392233])

#zernike_poly = Zernike_Combination((x,y), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
zernike_poly = Zernike_Combination((x,y), *coeffs)
zr = zernike_poly*pi
Kino = (Kino - zr)

Kino_img = ((Kino-np.min(Kino))/(np.max(Kino)-np.min(Kino)))*max_val

""" Tile-ing """
Nx = 1280
Ny = 768
x_reps = int(np.ceil(float(Nx)/float(sx)))+2
y_reps = int(np.ceil(float(Ny)/float(sy)))+2

kino_Final = np.zeros((sy*y_reps,sx*(x_reps)))

for i in xrange(0,x_reps):
    for j in xrange(0,y_reps):
        kino_Final[j*sy:(j+1)*sy,i*sx:(i+1)*sx] = Kino_img
        
noshift = kino_Final[0:Ny,0:Nx]
#shift_x = (int(sx*x_reps/2))-1280/2
#shift_y = (int(sy*y_reps - 768))/2
#kino_Final = kino_Final[shift_y:768+shift_y,shift_x:1280+shift_x]

#noshift = ((noshift-np.min(noshift))/(np.max(noshift)-np.min(noshift)))*255.0
phaseIm = Image.fromarray(noshift.astype(np.uint8))
