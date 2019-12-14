import numpy as np
import math
import matplotlib.pyplot as plt

epsilon = 3.084 * 10**(-3) #eV
sigma = 2.782 #Angstroms
sigma6 = sigma**6
sigma12 = sigma**12



def dotProduct(v1,v2):
    '''
    Arguments:

    v1,v2 - NumPy arrays of same size

    Returns:

    dp - result of dot product of v1 & v2
    '''

    dp = [x*y for x,y in zip(v1,v2)]
    dp = sum(dp)

    return dp


def LJ_potential(pos_vec):

    sig_dist = sigma/math.sqrt(dotProduct(pos_vec,pos_vec))

    return 4*epsilon * (sig_dist**12 - sig_dist**6)



def LJ_derivative(pos_vec):

    squareDist = dotProduct(pos_vec,pos_vec)

    term1 = 2 * sigma12 / squareDist**7
    term2 = sigma6 / squareDist**4


    return -24*epsilon * (term1 - term2) * pos_vec




def aGet(pos_vec,fprime,g):

    smolNum = 0.0001
    numerator = dotProduct(fprime(pos_vec),g)
    denomenator = dotProduct(fprime(pos_vec + smolNum * g)-fprime(pos_vec),g)
    return -smolNum * numerator / denomenator


def psuedoPBC(x1,x2):

    x1p, x2p = np.array([0.0,0.0,0.0]),np.array([0.0,0.0,0.0])

    #Make new coordinates same as old
    x1p =+ x1
    x2p += x2

    #Used to determine cutoff for pbc application
    x_size = 4



    if x1[0] - x2[0] > x_size/2 :
        x1p -= np.array([4.0,0.0,0.0])
    elif x1[0] - x2[0] < -x_size/2 :
        x1p += np.array([4.0,0.0,0.0])

    dif = x1p-x2p

    return dif



def getTotalPotential(distanceMatrix):
        '''
        Iterates through the distance matrix adding each pair of atoms contribution to the total potential energy of the lattice. This value is saved in the class instances attribute: latticePotential
        '''

        latticePotential = 0

        for row in distanceMatrix:
            latticePotential += LJ_potential(row[3])


def getDistanceMatrix(x):
        '''
        '''

        distanceMatrix = []

        atomCount, __  = x.shape
        for i in range(0,atomCount-1):
            for j in range(i+1,atomCount):

                #As this is the non pbc version, just take difference
                pos_vec = x[j,:] - x[i,:]

                distance = math.sqrt(dotProduct(pos_vec,pos_vec))


                if distance  <= self.cutoff:
                    distanceMatrix.append((i,j,pos_vec,distance))




def steepestDescent(g, f, x, a=0.1, h = 10**(-8)):

    #Get number of atoms for loops
    atomCount, __ =  x.shape

    

    #Initialise b_array to have all zeros so b's can be added
    b_array = np.zeros((x.shape))

    #Get distance matrix for total V and pos_vec's
    distMat = getDistanceMatrix(x)
    g0 = getTotalPotential(distMat)



    for row in distMat:

        #Unpacking atom indices and pos_vec from the distance matrix
        j,k,pos_vec = row[0:3]
        
        #gradient of current positions
        delg = f(pos_vec)
        
        #gauge of magnitude of gradients
        dg = math.sqrt(dotProduct(delg,delg))
        
        #scaling a for (i think) optimisation purposes
        b = a/dg 
        
        #Add to gradients np.array for addition to atoms array
        b_array[k] -= 0.5 * b*delg
        b_array[j] += 0.5 * b*delg




    count = 0
    for i in range(0,1):
        #Update all the atom's positions to reduce potential (hopefully)
        x += b_array
        for j in range(0,atomCount-1):
            for k in range(j,atomCount):


                #position vector of 1 w.r.t. 0
                pos_vec = x[k,:] - x[j,:]

                #Initial potential energy

                pot_dict{str(j)+str(k)} = g(pos_vec)

                #gradient of current positions
                delg = f(pos_vec)


                #gauge of magnitude of gradients
                dg = math.sqrt(dotProduct(delg,delg))

                #scaling a for (i think) optimisation purposes
                b = a/dg

                #Add to gradients np.array for addition to atoms array
                b_array[j] -= 0.5 * b*delg
                b_array[k] += 0.5 * b*delg

                #Calculate current potential to see if step size is too large or not
                g1 = g(pos_vec)
                if g1 > pot_dict[str(j)+str(k)]:
                    #If potential actually increased then step size is too large
                    a /= 2
                else:
                    #Update minimum potential if it has decreased
                    g0 = g1


       

        #position vector of 1 w.r.t. 0
        pos_vec = x[1,:] - x[0,:]

        #gradient of current positions
        delg = f(pos_vec)

        #gauge of magnitude of gradients
        dg =  math.sqrt(dotProduct(delg,delg))
        #scaling a for (i think) optimisation purposes
        b = a/dg

        #Calculate current potential to see if step size is too large or not
        g1 = g(pos_vec)
        if g1 > g0:
            #If potential actually increased then step size is too large
            a /= 2
            b = a/dg
        else:
            #Update minimum potential if it has decreased
            g0 = g1

        count += 1

    return pos_vec


x = np.array([[0.0,0.0,0.0],[1.0,0.0,0.0]])



result =steepestDescent(LJ_potential,LJ_derivative,x)


print("EOP")
