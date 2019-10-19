#Lennard-Jones potential for Ne crystal

import numpy as np
import math
import matplotlib.pyplot as plt

#Hard coded parameters for crystals
dimensionOfLattice = (2,1,1) #Equivalent to coding n's
LatticeConstant = 1

#Lennard Jones parameters

epsilon = 3.084 * 10**(-3) #ev
sigma = 2.782 #Angstroms





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
        self.cutoff = a + 0.001


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
        volume = np.dot(a1,np.cross(a2,a3))

        #Calcutlating the reciprocal vectors using the volume
        b1 = np.cross(a2,a3)/volume
        b2 = np.cross(a3,a1)/volume
        b3 = np.cross(a1,a2)/volume


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
        n1 = np.dot(b1,t)%1
        n2 = np.dot(b2,t)%1
        n3 = np.dot(b3,t)%1


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





    def getDistanceDict(self):
        '''
            NEED NEW EXPLANATION

        '''

        self.distanceCNT = {}


        atomCount, __  = self.atoms.shape

        for i in range(0,atomCount-1):
            for j in range(i+1,atomCount):

                #Apply PBC
                fracCord, PBCcoord = self.PBC(self.atoms[i,:],self.atoms[j,:])

                distance = np.linalg.norm(PBCcoord)

                if distance  <= self.cutoff:
                    found = False
                    for dist in sorted(self.distanceCNT.keys()):
                        if distance < dist +0.01 and distance > dist - 0.01:
                            self.distanceCNT[dist] += 1
                            found = True
                            break
                    if found == False:
                        self.distanceCNT[distance] = 0














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

        self.cutoff = (maxLength/1) * a + 0.01


        self.getDistanceDict()

        self.total_potential()


    def LJ_potential(self, distance):

        sig_dist = sigma/distance

        return 4*epsilon * (sig_dist**12 - sig_dist**6)


    def total_potential(self):

        forSum = [self.LJ_potential(x) * y for x,y in zip(self.distanceCNT.keys(),self.distanceCNT)]

        self.totalV = sum(forSum)









Ne1 = fcc("Ne",(2,2,2),LatticeConstant)
Ne2 = fcc("Ne",(3,3,3),LatticeConstant)
Ne3 = fcc("Ne",(4,4,4),LatticeConstant)
Ne4 = fcc("Ne",(5,5,5),LatticeConstant)




print(Ne1.totalV, Ne2.totalV,Ne3.totalV,Ne4.totalV)
