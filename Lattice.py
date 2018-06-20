import numpy as np
from PIL import Image

numX = 7
numY = 7
sX = 600
sY = 792

XGrid = np.linspace(0,sX,sX)
YGrid = np.linspace(0,sY,sY)
#
#xCenters = np.linspace(0+numX,sX-numX,numX)
#yCenters = np.linspace(0+numY,sY-numY,numY)
# 
#mu = 10
#
#
#
#intDis = np.zeros((sY,sX))
#
#for i in xrange(len(xCenters)):
#    for j in xrange(len(yCenters)):
#        intDis = intDis + np.exp( -((XGrid-xCenters[i])**2.0)[np.newaxis,:]/(2*mu**2)-((YGrid-yCenters[j])**2)[:,np.newaxis]/(2*mu**2) )
#        
#        
#latImg = (((intDis-np.min(intDis))/(np.max(intDis)-np.min(intDis)))*255).astype(np.uint8)
#latImg = Image.fromarray(latImg)    

points = np.array([4,5,6,7,6,5,4])
distance = 90
size = 85




intDis = np.zeros((sY,sX))

for i in range(len(points)):
    k = (len(points)-1)/2 -i
    print(k)
    yPos = sY/2 - distance*k
    for j in range(points[i]):
        u = (points[i]-1)/2 - j
        if(points[i]%2 == 0):
            xPos = sX/2 - distance*u - distance/2
        else:
            xPos = sX/2 - distance*u
        for ii in range(len(XGrid)):
            for jj in range(len(YGrid)):
                if ((XGrid[ii] - xPos)**2+(YGrid[jj] - yPos)**2 <= (size/2)**2):
                    intDis[jj,ii] = 1
                    
latImg = (((intDis-np.min(intDis))/(np.max(intDis)-np.min(intDis)))*255).astype(np.uint8)
latImg = Image.fromarray(latImg)