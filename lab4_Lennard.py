#Lennard-Jones potential for Ne crystal

import numpy as np
import math
import matplotlib.pyplot as plt
from collections import defaultdict
from functools import partial
import random

#Hard coded parameters for crystals
dimensionOfLattice = (2,2,2) #Equivalent to coding n's
LatticeConstant = 4.2

#Lennard Jones parameters

epsilon = 3.084 * 10**(-3) #eV
sigma = 2.782 #Angstroms
sigma6 = sigma**6
sigma12 = sigma**12


def writeToxyz(filename,crystal):
    ''' Takes some crystal class instance and writes its data to a .xyz file

    Arguments:

    filename -- Desired filename, should end in .xyz

    crystal -- instance of any crystal type (i.e. sc class or subclass of sc)


    '''
    crystalFile = open(filename,"wt")


    atomCountFinal = crystal.atoms.shape[0]


    crystalFile.write("{0}\n".format(atomCountFinal) )
    crystalFile.write("This is a {0} {1} \n".format(dimensionOfLattice,crystal.structure))
    for i in range(0,atomCountFinal):
        crystalFile.write("{0:s} \t {1:0.9f} \t {2:0.9f} \t {3:0.9f} \n".format(crystal.element,crystal.atoms[i,0],crystal.atoms[i,1],crystal.atoms[i,2]))

    return





def lineMinimisation(g,f_displaced,smolNum):

    numerator = np.sum(np.multiply(-g,g))
    denomenator = np.sum(np.multiply(f_displaced + g,g))
    return -smolNum * numerator / denomenator


def LJ_derivative(dist_PBC,pos_vec):

    
    term1 = 2 * sigma12 / dist_PBC**14 
    term2 = sigma6 / dist_PBC**8 


    return 24*epsilon * (term1 - term2) * pos_vec






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


def crossProduct(v1,v2):
    '''
    Arguments:

    v1,v2 - NumPy arrays of shape (1,3)

    Returns:

    v3 - Vector resulting from cross product
    '''


    pos1 = v1[1]*v2[2] - v1[2]*v2[1]
    pos2 = v1[2]*v2[0] - v1[0]*v2[2]
    pos3 = v1[0]*v2[1] - v1[1]*v2[0]

    v3 = np.array([pos1,pos2,pos3])


    return v3



class sc():
    '''Creates an instance of a simple cubic crystal

        Public methods:

        extendUnitCell() -- takes in an atom and returns an np.array containing
        atoms obtained by linear combination of the basis vectors within the
        dimension of the lattice


        Instance variables:

        self.element -- The chemical symbol of the element crystal atoms are

        self.structure -- The type of crystal structure, will change depending
        on subclass

        self.atoms -- An np.array containing all the atoms in the crystal. Each row is an atom while the columns are x,y and z postion respectively


        Constructor arguments:

        element -- The chemical symbol of the element crystal atoms are

        dimensionOfLattice -- A tuple containing the dimensions of the lattice in terms of unit cells

        a -- Lattice constant

        lVectors -- Lattice vectors for the given periodicity of unit cell

        recipVectors -- Reciprocal vectors to the lattice vectors

        cutoff -- Distance cutoff for nearest neighbour

        nearestN -- a list containing tuples that have
                    (index of 1st atom in self.atoms, index of 2nd atom, distance(1,2))
                    where 1 and 2 refer to lattice points.

    '''

    def __init__(self,element,dimensionOfLattice,a):

        self.element = element
        self.structure = "Simple Cubic"
        self.atoms = np.zeros((1,3))
        self.cutoff =  a + 0.001


        extraAtoms = self.extendUnitCell(self.atoms,dimensionOfLattice,a)
        self.atoms = extraAtoms


        #NEW CODE
        #Hardcoding the lattice vectors and the using a funciton
        #to find reciprocal vectors
        nx,ny,nz = dimensionOfLattice
        a1 = np.array((a*nx,0,0))
        a2 = np.array((0,a*ny,0))
        a3 = np.array((0,0,a*nz))
        self.lVectors = (a1,a2,a3)
        self.recipVectors, __ = self.getReciprocal()



    def extendUnitCell(self,atoms,dimensionOfLattice,a):
        ''' Returns all other similar atoms within the dimensions of the lattice


            This method will take in a np.array of atoms and return those and  all other atoms with the same fractional coordinates which are within the dimensions of the lattice



            Arguments:

            atoms -- an np.array of atoms (see sc class docstring) to be extended

            dimensionOfLattice -- A tuple containing the dimensions of the lattice  in terms of unit cells

            a -- Lattice constant
        '''

        unitVectors = [np.array([[a,0,0]]),np.array([[0,a,0]]),np.array([[0,0,a]])]



        for length, unitVector in zip(dimensionOfLattice,unitVectors):
            atomCount, __  = atoms.shape

            for i in range(0,length-1):
                #each time tempAtoms will be reading the atoms added in the
                #previous iteration
                tempAtoms = atoms[i*atomCount:(i+1)*atomCount,:] + unitVector
                atoms = np.concatenate((atoms,tempAtoms),axis = 0)


        return atoms



    def getReciprocal(self):
        '''
        '''
        #unpack lattice vectors
        a1,a2,a3 = self.lVectors

        #Finding the volume as directed in the lecture notes
        volume = dotProduct(a1,crossProduct(a2,a3))

        #Calcutlating the reciprocal vectors using the volume
        b1 = crossProduct(a2,a3)/volume
        b2 = crossProduct(a3,a1)/volume
        b3 = crossProduct(a1,a2)/volume


        #As requested in the lab notes returning both the reciprocal
        #vectors and the volume
        return (b1,b2,b3),volume





    #Make testing and non testing versions
    def PBC(self, l1,l2):
        '''
        takes in two lattice pointss and applies PBC where theorigin
        is the coordiante of the first lattice point

        will return fractional coordinates and the PBCcoordinates for
        the second lattice point (this is cartesian but in aproper)

        '''
        #get a's and b's

        a1,a2,a3 = self.lVectors
        b1,b2,b3 = self.recipVectors


        #Centering around l1
        t = l2-l1

        #Calculating fractional coordinates in line with lecturenotes
        n1 = dotProduct(b1,t)%1
        n2 = dotProduct(b2,t)%1
        n3 = dotProduct(b3,t)%1


        #A series of if statements to apply PBC i.e. move atoms
        #past the halfway point to equivalent location with PBC
        if n1 >= 0.5:
            n1 -= 1
        elif n1 < -0.5:
            n1 += 1

        if n2 >= 0.5:
            n2 -= 1
        elif n2 < -0.5:
            n2 += 1

        if n3 >= 0.5:
            n3 -= 1
        elif n3 < -0.5:
            n3  += 1


        #returnig the fracitnal coordinates and t2 as requested

        return (n1,n2,n3), n1*a1 + n2*a2 + n3*a3













class fcc(sc):
    '''Creates an instance of a face centred cubic crystal

        Subclassed from sc and extends that class by adding the extra atoms in
        a bcc unit cell

    '''

    def __init__(self,element,dimensionOfLattice,a):
        super().__init__(element,dimensionOfLattice,a)
        self.structure = "Face Centred Cubic"

        #Hard coded fcc atoms
        extraAtoms =np.array( [[0,0.5*a,0.5*a]+[0.5*a,0,0.5*a]+[0.5*a,0.5*a,0]])
        extraAtoms = extraAtoms.reshape((-1,3))
        extraAtoms = super().extendUnitCell(extraAtoms,dimensionOfLattice,a)
        self.atoms = np.concatenate((self.atoms,extraAtoms),axis = 0)


        maxLength = max(dimensionOfLattice)

        self.cutoff = a/math.sqrt(2) +0.001


    def LJ_potential(self, distance):

        sig_dist = sigma/distance

        return 4*epsilon * (sig_dist**12 - sig_dist**6)


    def getDistanceMatrix(self):
        '''
        Apply PBC to get a distance matrix containing all pairs of atoms closer
        than the cutoff.

        Creates:

        distanceMatrix - A list containing [atom 1 index, atom 2 index, distance, param string for paramDict]

        '''

        self.distanceMatrix = []

        atomCount, __  = self.atoms.shape
        #Changed to count every atom
        for i in range(0,atomCount):
            for j in range(0,atomCount):
                #Can optimise by only doing once through and then appending the
                #a list comprehension with elements 1,2 reversed
                if j != i:

                    #Apply PBC
                    fracCord, PBCcoord = self.PBC(self.atoms[i,:3],self.atoms[j,:3])

                    PBC_distance = math.sqrt(dotProduct(PBCcoord,PBCcoord))

                    if PBC_distance  <= self.cutoff:                
                        non_PBCcoord = self.atoms[j,:3] - self.atoms[i,:3]
                        norm_distance = math.sqrt(dotProduct(non_PBCcoord,non_PBCcoord))
                        self.distanceMatrix.append((i,j,PBC_distance,PBCcoord))


    def getTotalPotential(self):
        '''
        Iterates through the distance matrix adding each pair of atoms contribution to the total potential energy of the lattice. This value is saved in the class instances attribute: latticePotential
        '''

        self.latticePotential = 0

        for row in self.distanceMatrix:
            self.latticePotential += self.LJ_potential(row[2])





    def steepestDescent(self, fprime, h = 10**(-3), smolNum= 1):
        '''
        To avoid unnecessary calculations probably gonna remove element from distancematrix list
        '''
        for i in range(0,100):
            #Givesd correct shape for gradient
            g = np.zeros((self.atoms.shape))

            self.getDistanceMatrix()

            #Get negative of gradient
            for row in self.distanceMatrix:
                result = LJ_derivative(row[2],row[3]) / 2
                #Negative of gradient applied here
                g[row[0]] -= result
                g[row[1]] += result



            
            #Copy atoms before moving them for line minimisation
            temporary_copy = np.copy(self.atoms)
            self.atoms += sigma * g
            self.getDistanceMatrix()

            #Get f'(x + sigma*g)
            f_displaced = np.zeros((self.atoms.shape))
            for row in self.distanceMatrix:
                result = LJ_derivative(row[2],row[3]) / 2
                f_displaced[row[0]] += result
                f_displaced[row[1]] -= result

            alpha =  lineMinimisation(g,f_displaced,smolNum)

            print(alpha)
            self.atoms = temporary_copy + alpha * g
            




Ne = fcc("Ne",dimensionOfLattice,LatticeConstant)



#print(Ne.distanceMatrix)




writeToxyz("NeBeforeBef.xyz",Ne)

#Apply random perturbations to atoms
for i in range(0,150):
    RNG_atom = random.randint(0,Ne.atoms.shape[0]-1)
    RNG_coord = random.randint(0,2)
    temp = random.random() * 0.25
    Ne.atoms[RNG_atom,RNG_coord] += temp


writeToxyz("NeBefore.xyz",Ne)
Ne.steepestDescent(LJ_derivative)
writeToxyz("NeAfter.xyz",Ne)
