#Coulomb-Buckingham potential for KF crystal

import numpy as np
import math
import matplotlib.pyplot as plt

#Hard coded parameters for crystals
dimensionOfLattice = (5,5,5) #Equivalent to coding n's
LatticeConstant = 4.2 #Measured in Angstroms

#Coulomb Parameters
#In SI units
epsilon_not = 8.854 * 10**(-12)
e_charge = 1.602 * 10**(-19)

#Vaccum perittivity in eV * A * e^(-2)
ke = 1/(math.pi * 4 * epsilon_not) * e_charge * 10**10
#Effecive charge used
q = 1.7**2



#Buckingham parameters for MgO
#tuple contatins (A,rho,C) for anion <-> anion and cation <-> anion
#interactions
paramDict = {"mm":(4870.0,0.2670,77.0),"pm":(926.69,0.29909,0)}




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

        lVectors -- Lattice vectors for the given periodicity of unit cell

        recipVectors -- Reciprocal vectors to the lattice vectors

        cutoff -- Distance cutoff for nearest neighbour

        nearestN -- a list containing tuples that have
                    (index of 1st atom in self.atoms, index of 2nd atom, distance(1,2))
                    where 1 and 2 refer to lattice points


        Constructor arguments:

        element -- The chemical symbol of the element crystal atoms are

        dimensionOfLattice -- A tuple containing the dimensions of the lattice in terms of unit cells

        a -- Lattice constant

        init_pos -- Position of base atom which will be extended to
        desired peiodicity



    '''

    def __init__(self,element,dimensionOfLattice,a,init_pos):

        self.element = element
        self.structure = "Simple Cubic"
        self.atoms = np.array(init_pos)
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
        Calculates the reciprocal vectors for the lattice vectors
        for the particular class instance
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
        a fcc unit cell

    '''

    def __init__(self,element,dimensionOfLattice,a,init_pos):
        super().__init__(element,dimensionOfLattice,a,init_pos)
        self.structure = "Face Centred Cubic"

        init_pos = np.concatenate((init_pos,init_pos,init_pos),axis =1)
        #Hard coded fcc atoms
        extraAtoms =np.array( [[0,0.5*a,0.5*a]+[0.5*a,0,0.5*a]+[0.5*a,0.5*a,0]])
        extraAtoms += init_pos

        extraAtoms = extraAtoms.reshape((-1,3))
        extraAtoms = super().extendUnitCell(extraAtoms,dimensionOfLattice,a)
        self.atoms = np.concatenate((self.atoms,extraAtoms),axis = 0)

        #Add charge as 4th column
        atomCount, __ = self.atoms.shape

        if self.element == "Mg":
            chargeCol = np.zeros((atomCount,1)) + 2
        elif self.element == "O":
            chargeCol = np.zeros((atomCount,1)) -2

        self.atoms = np.concatenate((self.atoms,chargeCol),axis = 1)


        maxLength = max(dimensionOfLattice)

        self.cutoff = a * maxLength




    def Buckingham_potential(self,r,params):
        '''
        Calculates Buckingham potential of two atoms/ions

        Arguments:

        r - Seperation of the two atoms in Angstrom

        params - String containing "pm" if the potential is for cation <->
        anion and containgin "mm" if anion <-> anion

        Returns: Potential energy in eV
        '''

        A, rho, C = paramDict[params]

        return A * math.exp(-r/rho) - C/(r**6)


    def Coulomb_potential(self,r,params):
        '''
        Calculates Coulomb potential of two atoms/ions

        Arguments:

        r - Seperation of the two atoms in Angstrom

        params - String containing "pm" if the potential is for cation <->
        anion and containgin "mm" if anion <-> anion

        Returns: Potential energy in eV
        '''

        if params == "pm":
            q_product = - q
        else:
            q_product = q


        return (ke*q_product)/r


    def __add__(self,otherCrystal):
        '''
        Addition operator acting on an fcc crystal class instance now
        concatenates the array holding the atoms

        Has no return value - only a side effect
        '''

        self.atoms = np.concatenate((self.atoms,otherCrystal.atoms),axis = 0)

        return None


    def getDistanceMatrix(self):
        '''
        Apply PBC to get a distance matrix containing all pairs of atoms closer
        than the cutoff.

        Creates:

        distanceMatrix - A list containing [atom 1 index, atom 2 index, distance, param string for paramDict]

        '''

        self.distanceMatrix = []

        atomCount, __  = self.atoms.shape
        for i in range(0,atomCount-1):
            for j in range(i+1,atomCount):

                #Apply PBC
                fracCord, PBCcoord = self.PBC(self.atoms[i,:3],self.atoms[j,:3])

                distance = math.sqrt(dotProduct(PBCcoord,PBCcoord))


                if distance  <= self.cutoff:
                    if self.atoms[i,3] == self.atoms[j,3]:
                        if self.atoms[i,3] <0:
                            parameters = "mm"
                        else:
                            parameters = "pp"

                    else:
                        parameters = "pm"

                    self.distanceMatrix.append((i,j,distance,parameters))


    def getTotalPotential(self):
        '''
        Iterates through the distance matrix adding each pair of atoms contribution to the total potential energy of the lattice. This value is saved in the class instances attribute: latticePotential
        '''

        self.latticePotential = 0

        for row in self.distanceMatrix:

            #Get Buckingham potential
            if row[3] != "pp":
                self.latticePotential += self.Buckingham_potential(row[2],row[3])


            #Get Coulomn potential
            self.latticePotential += self.Coulomb_potential(row[2],row[3])



#Uncomment the below code to see the graph showing the lattice constant for
#which potential is mnimised
'''
lConstants = np.arange(4,4.5,0.05)
pots = []

for i in lConstants:
        Mgfcc = fcc("Mg",dimensionOfLattice,i,[[0,0,0]])
        Ofcc = fcc("O",dimensionOfLattice,i,[[0,0,i/2]])

        Mgfcc + Ofcc
        MgOfcc = Mgfcc

        MgOfcc.getDistanceMatrix()
        MgOfcc.getTotalPotential()

        pots.append(MgOfcc.latticePotential)



print("The below graph shows the minimum occuring at 4.2 Angstroms as expected for MgO")
plt.plot(lConstants,pots)
plt.ylabel("Potential (eV)")
plt.xlabel("Lattice Constant (Angstrom)")
plt.title("Potential Minimising Lattice Constant")
plt.show()
'''





#To see test for converging average energy per per atom uncomment belwo section
'''
dims = [(x,x,x) for x in range(1,6)]
pots = []


for i in dims:
        Mgfcc = fcc("Mg",i,LatticeConstant,[[0,0,0]])
        Ofcc = fcc("O",i,LatticeConstant,[[0,0,LatticeConstant/2]])

        Mgfcc + Ofcc
        MgOfcc = Mgfcc

        MgOfcc.getDistanceMatrix()
        MgOfcc.getTotalPotential()

        pots.append(MgOfcc.latticePotential/(i[1]**3))

print(Below is the average energy per atom for cubes of size 1x1x1 to 5x5x5: \n)
print(pots)
'''

Mgfcc = fcc("Mg",dimensionOfLattice,LatticeConstant,[[0,0,0]])
Ofcc = fcc("O",dimensionOfLattice,LatticeConstant,[[0,0,LatticeConstant/2]])

Mgfcc + Ofcc
MgOfcc = Mgfcc

MgOfcc.getDistanceMatrix()
MgOfcc.getTotalPotential()

print("Given an MgO lattice of size {0} and with lattice constant {1} Angstroms the total potential was found to be {2} eV".format(dimensionOfLattice,LatticeConstant,MgOfcc.latticePotential))
