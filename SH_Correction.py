# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 11:36:34 2016

@author: thomas
"""

import numpy as np
import scipy.optimize as opt
from scipy.misc import derivative as der
import matplotlib.pyplot as plt
from PIL import Image
from skimage.draw import line
import time 

def Zernike_Combination(x,y, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15):
    rho = np.zeros((np.size(y),np.size(x)))    
    phi = np.zeros((np.size(y),np.size(x)))    
    for i in range(0,np.size(x)):
        for j in range(0,np.size(y)):
            rho[j,i] = (x[i]**2.0+y[j]**2.0)**0.5
            phi[j,i] = np.arctan2(y[j],x[i])
    p1 = 1.0 #Piston
    
    
    p2 = 2.0*rho*np.cos(phi)#Tip (x-tilt)
    
   
    p3 = 2.0*rho*np.sin(phi) #Tilt (y-tilt)
    
    
    p4 = (3.0**0.5)*(2.0*rho**2.0 - 1.0)#Defocus
    
   
    p5 = (6.0**0.5)*(rho**2.0)*np.sin(2.0*phi) #Oblique Astigmatism
    
     
    p6 = (6.0**0.5)*(rho**2.0)*np.cos(2.0*phi) #Vertical Astigmatism  
    
    
    p7 = (8.0**0.5)*(3.0*rho**3.0 - 2.0*rho)*np.sin(phi) #Vertical Coma   
    
       
    p8 = (8.0**0.5)*(3.0*rho**3.0 - 2.0*rho)*np.cos(phi)#Horizontal Coma 
    
       
    p9 = (8.0**0.5)*(rho**3.0)*np.sin(3.0*phi)#Vertical Trefoil 
    
    
    p10 = (8.0**0.5)*(rho**3.0)*np.cos(3.0*phi) #Oblique Trefoil
    
    #Primary Spherical    
    p11 = (5.0**0.5)*(6.0*rho**4.0 - 6.0*rho**2.0 + 1.0)
    
    #Vertical Secondary Astigmatism    
    p12 = (10.0**0.5)*(4*rho**4 - 3*rho**2)*np.cos(2*phi)
    
    #Horizontal Secondary Asstigmatism
    p13 = (10.0**0.5)*(4.0*rho**4.0 - 3.0*rho**2.0)*np.sin(2.0*phi)
    
    #Vertical Quadrafoil    
    p14 = (10.0**0.5)*(rho**4.0)*np.cos(4.0*phi)
    
    #Oblique Quadrafoil
    p15 = (10.0**0.5)*(rho**4.0)*np.sin(4.0*phi)
    
    full = m1*p1 + m2*p2 + m3*p3 + m4*p4 + m5*p5 + m6*p6 + m7*p7 + m8*p8 + m9*p9 + m10*p10 + m11*p11 + m12*p12 + m13*p13 + m14*p14 + m15*p15    
    
    return full.ravel()

def Zernike_Noll(j, x,y):
    rho = np.zeros((np.size(y),np.size(x)))    
    phi = np.zeros((np.size(y),np.size(x)))    
    for i in range(0,np.size(x)):
        for n in range(0,np.size(y)):
            rho[n,i] = (x[i]**2.0+y[n]**2.0)**0.5
            phi[n,i] = np.arctan2(y[n],x[i])
    p = np.zeros((np.size(y), np.size(x)))
    print(j)    
    if (j == 1):    
        #Piston
        p = np.ones((np.size(y), np.size(x))) 
    elif (j == 2):
        #Tip (x-tilt)
        p = 2.0*rho*np.cos(phi)
    elif (j == 3):
        #Tilt (y-tilt)
        p = 2.0*rho*np.sin(phi)
    elif (j == 4):
        #Defocus
        p = (3.0**0.5)*(2.0*rho**2.0 - 1.0)
    elif (j == 5):
        #Oblique Astigmatism
        p = (6.0**0.5)*(rho**2.0)*np.sin(2.0*phi)
    elif (j == 6):
        #Vertical Astigmatism    
        p = (6.0**0.5)*(rho**2.0)*np.cos(2.0*phi)
    elif (j == 7):
        #Vertical Coma    
        p = (8.0**0.5)*(3.0*rho**3.0 - 2.0*rho)*np.sin(phi)
    elif (j == 8):
        #Horizontal Coma    
        p = (8.0**0.5)*(3.0*rho**3.0 - 2.0*rho)*np.cos(phi)
    elif (j == 9):
        #Vertical Trefoil    
        p = (8.0**0.5)*(rho**3.0)*np.sin(3.0*phi)
    elif (j == 10):
        #Oblique Trefoil
        p = (8.0**0.5)*(rho**3.0)*np.cos(3.0*phi) 
    elif (j == 11):
        #Primary Spherical    
        p = (5.0**0.5)*(6.0*rho**4.0 - 6.0*rho**2.0 + 1.0)
    elif (j == 12):
        #Vertical Secondary Astigmatism    
        p = (10.0**0.5)*(4*rho**4 - 3*rho**2)*np.cos(2*phi)
    elif (j == 13):
        #Horizontal Secondary Asstigmatism
        p = (10.0**0.5)*(4.0*rho**4.0 - 3.0*rho**2.0)*np.sin(2.0*phi)
    elif (j == 14):
        #Vertical Quadrafoil    
        p = (10.0**0.5)*(rho**4.0)*np.cos(4.0*phi)
    elif (j == 15):
        #Oblique Quadrafoil
        p = (10.0**0.5)*(rho**4.0)*np.sin(4.0*phi)
    else:
        print("Noll index out of bounds")
    return p
 
def Displacement_Map(data_holo,holo_x,holo_y,lat_x,lat_y):
    f = 0.01
    sy,sx = data_holo.shape
    dispalcement_map  = np.zeros((sy,sx)) 
    gradients = np.zeros(2*len(lat_x))
    x = np.linspace(-1,1,sx)
    y = np.linspace(-1,1,sy)
    for i in range(len(holo_x)):
        rr, cc = line(int(holo_y[i]), int(holo_x[i]), int(lat_y[i]), int(lat_x[i]))
        val = (((holo_x[i]-lat_x[i])**2.0 + (holo_y[i]-lat_y[i])**2.0)**0.5)
        sign_val = np.sign((holo_y[i]-lat_y[i])/(holo_x[i]-lat_x[i]))
        #dispalcement_map[rr, cc] =  255
        dispalcement_map[rr, cc] = val*sign_val
        gradients[i] = (x[int(lat_x[i])] - x[int(holo_x[i])])/f
        gradients[i+len(lat_x)] = (y[int(lat_y[i])] - y[int(holo_y[i])])/f
    return dispalcement_map, gradients   
   
def Derrivative_Matrix(lat_x, lat_y, data_lat):
    l = 2*np.size(lat_x)
    z_der = np.zeros((l,15))
    (sy,sx) = np.shape(data_lat)    
    x = np.linspace(-1,1,sx)
    y = np.linspace(-1,1,sy)
    dx = x[2]-x[1]
    dy = y[2]-y[1]
    
    for i in range(15):
        noll = i+1
        zp = Zernike_Noll(noll,(x,y))
        xgrad = np.gradient(zp,dx)[1]
        ygrad = np.gradient(zp,dy)[0]
        for j in range(len(lat_x)):
            z_der[j,i] = xgrad[lat_y[j],lat_x[j]]
            z_der[j+len(lat_x),i] = ygrad[lat_y[j],lat_x[j]] 
    
    return z_der


def Derivative_Fit(z_der,m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15):
    coeffs = np.array([m1,m2,m3,m4,m5,m6,m7,m8,m9,m10, m11, m12, m13, m14, m15])
    grads = np.dot(z_der,coeffs)
    return grads
    

def Get_Max_Points(data,num):
    x = np.zeros(num)
    y = np.zeros(num)
    for i in range(num):
        figMan = plt.get_current_fig_manager()
        figMan.window.showMaximized()
        plt.imshow(data)
        points = plt.ginput(2)
        plt.close()
        points = np.floor(points)
        roi = data[points[0][1]:points[1][1],points[0][0]:points[1][0]]
        roi = roi/np.max(roi)
        if(np.max(roi) == 255):
            print('Oh no...')
        total = np.sum(np.sum(roi,1),0)
        index_array = np.linspace(0, np.shape(roi)[1]-1,np.shape(roi)[1])
        x_loc = (1/total)*np.sum(np.sum(roi,0)*index_array)        
        index_array = np.linspace(0, np.shape(roi)[0]-1,np.shape(roi)[0]) 
        y_loc = (1/total)*np.sum(np.sum(roi,1)*index_array)
#        print ('max_local: ', max_index_loc)        
#        max_index = (points[0][1]+max_index_loc[0],points[0][0]+max_index_loc[1])
#        print('max_index: ', max_index)
        x[i] = points[0][0]+x_loc
        y[i] = points[0][1]+y_loc
    return x,y    

def Get_Centers(data, col_numb, cut, length_cut):
    x = np.zeros(0)
    y = np.zeros(0)
        
    for i in range(col_numb):
        """ Choose Column"""        
        figMan = plt.get_current_fig_manager()
        figMan.window.showMaximized()
        plt.imshow(data)
        points = plt.ginput(2)
        plt.close()
        points = np.floor(points)
        roi = data[points[0][1]:points[1][1],points[0][0]:points[1][0]]
        """Create threshold array"""
        roi = np.sum(roi,1)
        roi = roi/(np.max(roi))
        roi = (roi>cut).astype(int)
        
        """Find spot locations"""
        target = 1   
        j = 0         
        up = np.zeros(0)
        down = np.zeros(0)        
        for i in range(np.size(roi)):
            if(roi[i]==target and target == 1):
                j = i            
                target = 0
            elif(roi[i]==target and target == 0):
                #down = np.append(down,i+points[0][1])
                if (i-j>length_cut):
                    up = np.append(up,j+points[0][1])
                    down = np.append(down,i+points[0][1])
                target = 1
        """Calculate spot center"""
        print(np.size(up), np.size(down))
        time.sleep(0.5)
        for i in range(np.size(up)):
            roi = data[up[i]:down[i], points[0][0]:points[1][0]]
            roi = roi/np.max(roi)
            if(np.max(roi) == 255):
                print('Oh no...')
            total = np.sum(np.sum(roi,1),0)
            index_array = np.linspace(0, np.shape(roi)[1]-1,np.shape(roi)[1])
            x_loc = (1.0/total)*np.sum(np.sum(roi,0)*index_array)        
            index_array = np.linspace(0, np.shape(roi)[0]-1,np.shape(roi)[0]) 
            y_loc = (1.0/total)*np.sum(np.sum(roi,1)*index_array)
            #print ('max_local: ', max_index_loc)        
            #max_index = (points[0][1]+max_index_loc[0],points[0][0]+max_index_loc[1])
            #print('max_index: ', max_index)
            x = np.append(x, points[0][0] + x_loc)
            y = np.append(y, y_loc + up[i])
    
    return x,y

def Get_Data(name_holo, name_lattice):
    im = Image.open(name_holo).convert('L')
    raw_data = np.array(im.getdata()).astype(float)
    sx,sy = im.size
    data_holo = np.reshape(raw_data,(sy,sx))

    im = Image.open(name_lattice).convert('L')
    raw_data = np.array(im.getdata()).astype(float)
    sx,sy = im.size
    data_lat = np.reshape(raw_data,(sy,sx))

    plt.imshow(data_lat)
    print("Choose Horizontal, then vertical.")
    points = plt.ginput(4)
    plt.close()
    np.floor(points)
    
    data_holo = data_holo[points[2][1]:points[3][1],points[0][0]:points[1][0]]
    data_lat = data_lat[points[2][1]:points[3][1],points[0][0]:points[1][0]]

    return data_holo, data_lat

def Remove_Zero_Order(data_holo, data_lat):
    
    figMan = plt.get_current_fig_manager()
    figMan.window.showMaximized()       
    plt.imshow(data_holo)
    points = plt.ginput(2)
    plt.close()
    np.floor(points)
    
    data_holo[points[0][1]:points[1][1],points[0][0]:points[1][0]] = 0

    figMan = plt.get_current_fig_manager()
    figMan.window.showMaximized()    
    plt.imshow(data_lat)
    points = plt.ginput(2)
    plt.close()
    np.floor(points)
    
    data_lat[points[0][1]:points[1][1],points[0][0]:points[1][0]] = 0

def Do_Zernike_Fit(disp, z_der, initial_guess, data_holo):
    sy, sx = data_holo.shape     
    
    x = np.linspace(-1,1,sx)
    y = np.linspace(-1,1,sy)
     
    #z_der = z_der.ravel()
    popt, pcov = opt.curve_fit(Derivative_Fit, z_der, disp, p0=initial_guess)

    zernike = Zernike_Combination((x,y),*popt)
    zr = zernike.reshape(sy,sx)
    
    nx = 600
    ny = 792
    x = np.linspace(-nx/2,nx/2,nx)/(sx/2)
    y = np.linspace(-ny/2,ny/2,ny)/(sy/2)    
    
    zernike = Zernike_Combination((x,y),*popt)
    zr_full = zernike.reshape(ny,nx)
    
    return zr, zr_full, popt, pcov

def Convert_To_Image(data):
    data = ((data-np.min(data))/(np.max(data)-np.min(data)))*255.0
    img = Image.fromarray(data.astype(np.uint8))
    return img
    
figMan = plt.get_current_fig_manager()
#figMan.window.showMaximized()

number_points = 68
number_cols = 11

data_holo, data_lat = Get_Data('holo_500.png','pin_hole_500.png')
print("Loaded Data")
#Remove_Zero_Order(data_holo, data_lat)

#x = np.zeros(0)
#y = np.zeros(0)
    
#for i in xrange(number_cols):
#    """ Choose Column"""        
#    figMan = plt.get_current_fig_manager()
#    figMan.window.showMaximized()
#    plt.imshow(data_lat)
#    points = plt.ginput(2)
#    plt.close()
#    points = np.floor(points)
#    roi = data_lat[points[0][1]:points[1][1],points[0][0]:points[1][0]]
#    """Create threshold array"""
#    rois = np.sum(roi,1)
#    rois = rois/np.max(rois)
#    rois = (rois>0.07).astype(int)
#    
#    """Find spot locations"""
#    target = 1   
#    j = 0         
#    up = np.zeros(0)
#    down = np.zeros(0)        
#    for i in xrange(np.size(rois)):
#        if(rois[i]==target and target == 1):
#            j = i            
#            target = 0
#        elif(rois[i]==target and target == 0):
#            #down = np.append(down,i+points[0][1])
#            if (i-j>7):
#                up = np.append(up,j+points[0][1])
#                down = np.append(down,i+points[0][1])
#            target = 1
#            
#    """Calculate spot center"""
#    print(np.size(up), np.size(down))
#    
#    for i in xrange(np.size(up)):
#        roi2 = data_lat[up[i]:down[i], points[0][0]:points[1][0]]
#        roi2 = roi2/np.max(roi2)
#        if(np.max(roi2) == 255):
#            print('Oh no...')
#        total = np.sum(np.sum(roi,1),0)
#        index_array = np.linspace(0, np.shape(roi2)[1]-1,np.shape(roi2)[1])
#        x_loc = (1.0/total)*np.sum(np.sum(roi2,0)*index_array)        
#        index_array = np.linspace(0, np.shape(roi2)[0]-1,np.shape(roi2)[0]) 
#        y_loc = (1.0/total)*np.sum(np.sum(roi2,1)*index_array)
#        #print ('max_local: ', max_index_loc)        
#        #max_index = (points[0][1]+max_index_loc[0],points[0][0]+max_index_loc[1])
#        #print('max_index: ', max_index)
#        x = np.append(x, points[0][0]+x_loc)
#        y = np.append(y, points[0][1]+y_loc)

print("Find centers for pin hole image")
center_lat_x, center_lat_y = Get_Centers(data_lat,number_cols, 0.07, 8)
print("Find centers for hologram image")
center_holo_x ,center_holo_y = Get_Centers(data_holo, number_cols, 0.1, 10)

#max_holo_x ,max_holo_y = Get_Max_Points(data_holo,number_points)
#max_lat_x, max_lat_y = Get_Max_Points(data_lat,number_points)

print("Center of points calculated, generating dispalcement map")

dis_map, disp = Displacement_Map(data_holo,center_holo_x, center_holo_y,center_lat_x,center_lat_y)
print("Generating derivative matrix around center points")
z_der = Derrivative_Matrix(center_lat_x, center_lat_y, data_lat)
print("Doing fit using Zernike polynomials")
initial_guess = (0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
zernike, zr_full, coeffs, pcov = Do_Zernike_Fit( disp, z_der, initial_guess, data_holo)
zernike_image = Convert_To_Image(zernike)
coeffs_abs = np.abs(coeffs)