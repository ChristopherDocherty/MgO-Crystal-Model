# In[]:

import numpy as np
import math

#Hard coded parameters for crystals
dimensionOfLattice = (1,1,1) #Equivalent to coding n's
a= 5

t = (16.5,12.3,17.1)


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

#Part 2 find fractinal coordinates
# In[]:

def getFracCoord(a1,a2,a3):


    #find vectors ni for calculation 
    n1 = np.dot(b1,t)%1 - 0.5
    n2 = np.dot(b2,t)%1 - 0.5
    n3 = np.dot(b3,t)%1 - 0.5

    
    #Sum the vectors ai with ni prefactors
    return n1*a1 + n2*a2 + n3*a3




# In[]:



scCutoff = a
fccCutoff = math.sqrt(2)/4 * a
bccCutoff = math.sqrt(3)/16 * a


#use np.linalg.norm() Just put in difference vector


#Second nicer way is to store in both each time and then cnosider one less atom each time n! checks v.s. n^n checks

#will use np matrix to do the same as in the powerpoint except label1 and label2 will be three x,y,z coordinates of the atoms
#in question HERE


#Use self.atoms for number and then iterate over the row index
#Again iterate over row index but one greater this time
#Check if np.linalg.norm() is less than cutoff
#If yes then copy in row index for each and put in distance
#After fully iterating just use list comprehension or somethign to fip the first two columns and add them in

#use np.flip to make a copy conttenate distances then concatenate that to the original matrix



#Check should give the correct number of neighbouring atoms

#COuld potentially save neighbours to xyz file and have that serve as proof of validity














#%%
