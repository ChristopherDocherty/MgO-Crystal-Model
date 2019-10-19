#Program to produce .xyz files for Sc, Bcc, Fcc, and
#Diamond crystalline structures with the intent to display
#in VMD

import numpy as np
import math


#Hard coded parameters for crystals
dimensionOfLattice = (4,4,1) #Equivalent to coding n's
LatticeConstant = 1



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

        #Create an empty list to store neighbour list and then call funciton  to populate it
        self.nearestN = []
        self.findNearest()


    def extendUnitCell(self,atoms,dimensionOfLattice,a):
        ''' Returns all other similar atoms within the dimensions of the lattice


            This method will take in a np.array of atoms and return those and  all other atoms with the same fractional coordinates which are within the dimensions of the lattice



        Arguments:

        atoms -- an np.array of atoms (see sc class docstring) to be extended

        dimensionOfLattice -- A tuple containing the dimensions of the lattice in terms of unit cells

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
        takes in two lattice pointss and applies PBC where the origin
        is the coordiante of the first lattice point

        will return fractional coordinates and the PBC coordinates for
        the second lattice point (this is cartesian but in a proper)

        '''
        #get a's and b's

        a1,a2,a3 = self.lVectors
        b1,b2,b3 = self.recipVectors


        #Centering around l1
        t = l2-l1

        #Calculating fractional coordinates in line with lecture notes
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





    def findNearest(self):
        '''
        Goes through each pair of atoms only once adn calculates
        teh distance with PBC, if it is suficiently low to be a nearest neighbour the two atoms' indexes in self.atoms
        are saved into self.nearestN with their corresponding distances

        This list will contain a copy of all the pairs of atoms
        but to be in line with the notes all permutations of
        aotms are saved not just all the combinations


        this ensures that each atom in row 1 will have the correct number of nearest neighbours in row 2 and that
        the lenght of the list is:
        #unit cells * coordination number * lattice points per unit cell

        '''

        atomCount, __  = self.atoms.shape
        for i in range(0,atomCount-1):
            for j in range(i+1,atomCount):

                #Apply PBC
                fracCord, PBCcoord = self.PBC(self.atoms[i,:],self.atoms[j,:])

                distance = np.linalg.norm(PBCcoord)

                if distance  <= self.cutoff:
                    self.nearestN.append((i,j,distance))






        #In order to save all permutations tall the elements
        #of the list must be copied with the two atoms reversed
        #This is done through a list comprehension with tuple
        #flipping
        forExtend = [(x[1],x[0],x[2]) for x in self.nearestN]

        #Adding the extra permutations onto the
        self.nearestN.extend(forExtend)
        #Sort according to the first row atom
        self.nearestN.sort(key = lambda a : a[0])





class bcc(sc):
    '''Creates an instance of a body centred cubic crystal

        Subclassed from sc and extends that class by adding the extra atoms in
        a bcc unit cell

    '''

    def __init__(self,element,dimensionOfLattice,a):
        super().__init__(element,dimensionOfLattice,a)
        self.structure = "Body Centred Cubic"


        #Hard coded bcc atom
        extraAtoms = 0.5*np.array([[a,a,a]])
        extraAtoms = super().extendUnitCell(extraAtoms,dimensionOfLattice,a)
        self.atoms = np.concatenate((self.atoms,extraAtoms),axis = 0)



        #Add to the true value to cunter floating point error
        self.cutoff = a*math.sqrt(3)/2 + 0.01





        self.nearestN = []
        self.findNearest()



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


        self.cutoff = a/math.sqrt(2) +0.001


        self.nearestN = []
        self.findNearest()






#For part 1 of lab 2:




a = LatticeConstant

print("the lattice constant is {} \n".format(a))


#Calculating the colume and reciprocal vectors for SC primitive
part1Crystal = sc("Si",(1,1,1),a)

#In the case of simple cubci the unit and primitive cell are identical hence
recVec, V = part1Crystal.getReciprocal()

print("For a primitive simple cubic cell with the chosen lattice constant, the reciprocal vectors are {0}, {1} and {2} and the volume is {3} \n".format(recVec[0],recVec[1],recVec[2],V))





#Calculating the volume and reciprocal vectors for bcc primitive

t1 = a/2 *np.array([-1,1,1])
t2 = a/2 *np.array([1,-1,1])
t3 = a/2 *np.array([1,1,-1])


part1Crystal.lVectors =(t1,t2,t3)

recVec, V = part1Crystal.getReciprocal()

print("For a primitive body centred cubic cell with the chosen lattice constant, the reciprocal vectors are {0}, {1} and {2} and the volume is {3} \n".format(recVec[0],recVec[1],recVec[2],V))





#Calculating the volume and reciprocal vectors for fcc primitive

t1 = a/2 *np.array([0,1,1])
t2 = a/2 *np.array([1,0,1])
t3 = a/2 *np.array([1,1,0])


part1Crystal.lVectors =(t1,t2,t3)

recVec, V = part1Crystal.getReciprocal()

print("For a primitive face centred cubic cell with the chosen lattice constant, the reciprocal vectors are {0}, {1} and {2} and the volume is {3} \n".format(recVec[0],recVec[1],recVec[2],V))



#Part 2 of the lab is the method PBC() for sc and associated subclasses

#If you need to test specific if PBC() works for specific instances
#then it can be called as a class method. lattice and reciprocal
#vectors can be specified as in part 1



#Part 3 of the lab is contained in findNearest() class method



#Part 4 of the lab is shown  below

sc1 = sc("Si",(4,4,1),LatticeConstant)

print("below is sc1's nearest neighbour list of expected size {} (1 lattice point * 16 unit cells * 4 NN)".format(len(sc1.nearestN)))
print(sc1.nearestN)
print("\n")



bcc1 = bcc("Si",(3,3,3),LatticeConstant)

print("below is bcc1's nearest neighbour list of expected size {} (2 lattice point * 27 unit cells * 8 NN - as this is a 2d lattice like the powerpoint)".format(len(bcc1.nearestN)))
print(bcc1.nearestN)
print("\n")


fcc1 = fcc("Si",(3,3,3),LatticeConstant)
print("below is fcc1's nearest neighbour list of expected size {} (4 lattice point * 27 unit cells *  12NN)".format(len(fcc1.nearestN)))
print(fcc1.nearestN)
print("\n")

#All the above parameters can be altered if a different test is required
