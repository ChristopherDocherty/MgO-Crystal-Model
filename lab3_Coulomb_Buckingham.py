#Coulomb-Buckingham potential for KF crystal

import numpy as np
import math
import matplotlib.pyplot as plt

#Hard coded parameters for crystals
dimensionOfLattice = (4,4,4) #Equivalent to coding n's
LatticeConstant = 4.2 #Measured in Angstroms

#Coulomb Parameters

epsilon_not = 8.854 * 10**(-12)
e_charge = 1.602 * 10**(-19)



#Buckingham parameters

#MgO
#F- <-> F- interactions
Amm = 4870.0 #eV
rhomm =  0.2670 # Angstroms
Cmm = 77.0 #Angstrom ^(-6)

#K+ <-> F- interactions
Apm = 929.69 #eV
rhopm = 0.29909  #Angstroms
Cpm = 0 #Angstrom ^(6) eV
'''

#KCl
#K+ <-> K+ interactions
App = 3796.9#3796.9 #eV
rhopp = 0.2603#0.2603 #Angstroms
Cpp = 52.0 #Angstrom ^(-6)

#F- <-> F- interactions
Amm = 1227.2 #eV
rhomm = 0.3214 # Angstroms
Cmm = 165.4 #Angstrom ^(-6)

#K+ <-> F- interactions
Apm = 4117.9 #eV
rhopm = 0.304  #Angstroms
Cpm = 124.9 #Angstrom ^(6) eV
'''



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
                        self.distanceCNT[distance] = 1














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


        maxLength = max(dimensionOfLattice)



        self.cutoff = a/2


        self.getDistanceDict()

        self.total_potential()


    def Buckingham_potential(self,r,two_atom):

        if two_atom == True:
            A, rho, C = (Apm,rhopm,Cpm)
        else:
            A, rho, C = (Amm,rhomm,Cmm)
        temp=A * math.exp(-r/rho)
        temp2= - C/(r**6)
        return temp+temp2


    def Coulomb_potential(self,r,two_atom):

        qplus = 2
        qminus = -2

        r = r


        if two_atom == True:
            q_product = qplus*qminus
        else:
            if self.element == "K":
                q_product = qplus**2
            elif self.element == "F":
                q_product = qminus**2

        #ke = 1/(math.pi * 4 * epsilon_not)
        ke = 14.3981
        return (ke*q_product)/r



    def total_potential(self):

        if self.element == "K":
            forSum = [self.Coulomb_potential(x,False) * y for x,y in zip(self.distanceCNT.keys(),self.distanceCNT)]
        else:
            forSum = [(self.Coulomb_potential(x,False) + self.Buckingham_potential(x,False)) * y for x,y in zip(self.distanceCNT.keys(),self.distanceCNT)]

        self.totalV = sum(forSum)


    def get_combined_DistanceDict(self,other):

        self.combineDistanceCNT = {}

        atomCount, __  = self.atoms.shape

        for i in range(0,atomCount):
            for j in range(0,atomCount):

                #Apply PBC
                fracCord, PBCcoord = self.PBC(self.atoms[i,:],other.atoms[j,:])


                distance = np.linalg.norm(PBCcoord)

                if distance  <= self.cutoff:
                    found = False
                    for dist in sorted(self.combineDistanceCNT.keys()):
                        if distance < 1.01* dist and distance > 0.99*dist:
                            self.combineDistanceCNT[dist] += 1
                            found = True
                            break
                    if found == False:
                        self.combineDistanceCNT[distance] = 1


    def two_atom_basis_potential(self,other):

        self.get_combined_DistanceDict(other)



        forSum = [(self.Coulomb_potential(x,True) +self.Buckingham_potential(x,True))*y for x,y in zip(self.combineDistanceCNT.keys(),self.combineDistanceCNT)]

        return sum(forSum)




both_potential = []
singles_potential = []

rs =  np.arange(4.0,10.0,0.3)
dims = [(x,x,x) for x in range(1,6)]


for x in rs:

    #Kfcc = fcc("K",x,LatticeConstant,np.zeros((1,3)))
    #Ffcc = fcc("F",x,LatticeConstant,))

    Kfcc = fcc("K",dimensionOfLattice,x,np.zeros((1,3)))
    Ffcc = fcc("F",dimensionOfLattice,x,np.array([[0,0,x/2]]))


    plus_minus_potential = Kfcc.two_atom_basis_potential(Ffcc)


    both_potential.append(plus_minus_potential + Kfcc.totalV + Ffcc.totalV)
    print(plus_minus_potential,Kfcc.totalV,Ffcc.totalV)

plt.plot(rs,both_potential)
plt.show()


#ePerAtom = [both_potential[x]/((x**3)*2) for x in range(2,6)]
#print(ePerAtom)
