# In[]:

import numpy as np
import math

#Hard coded parameters for crystals
dimensionOfLattice = (1,1,1) #Equivalent to coding n's
a= 2

t = np.array([[16.5,12.3,17.1],[16.6,12.4,17.2]])



nx,ny,nz = dimensionOfLattice

a1 = np.array((a*nx,0,0))
a2 = np.array((0,a*ny,0))
a3 = np.array((0,0,a*nz))


def getReciprocal(a1,a2,a3):

    b1 = np.cross(a2,a3)/(np.dot(a1,np.cross(a2,a3)))
    b2 = np.cross(a3,a1)/(np.dot(a1,np.cross(a2,a3)))
    b3= np.cross(a1,a2)/(np.dot(a1,np.cross(a2,a3)))

    return b1,b2,b3



b1,b2,b3 = getReciprocal(a1,a2,a3)


def PBC(l1,l2):
    '''
    takes in two lattice pointss and applies PBC where the origin
    is the coordiante of the first lattice point

    will return fractional coordinates and the PBC coordinates for
    the second lattice point (this is cartesian but in a proper)

    '''
    t = l2-l1


    n1 = np.dot(b1,t)%1
    n2 = np.dot(b2,t)%1
    n3 = np.dot(b3,t)%1


    if n1 > 0.5:
        n1 -= 1
    elif n1 < -0.5:
        n1 += 1

    if n2 > 0.5:
        n2 -= 1
    elif n2 <-0.5:
        n2 += 1

    if n3 > 0.5:
        n3 -= 1
    elif n3 < -0.5:
        n3  += 1


    return (n1,n2,n3), n1*a1 + n2*a2 + n3*a3





l1 = np.array([0.3,0.3,0.3])
l2 = np.array([15.5,16.3,14.2])

print(PBC(l1,l2))











#%%
