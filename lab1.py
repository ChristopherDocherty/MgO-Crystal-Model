#Program to produce .xyz files for Sc, Bcc, Fcc, and
#Diamond crystalline structures with the intent to display
#in VMD

import numpy as np


#Hard coded parameters for crystals
dimensionOfLattice = (2,2,2)
LatticeConstant = 5



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

    '''

    def __init__(self,element,dimensionOfLattice,a):

        self.element = element
        self.structure = "Simple Cubic"
        self.atoms = np.zeros((1,3))

        extraAtoms = self.extendUnitCell(self.atoms,dimensionOfLattice,a)

        self.atoms = extraAtoms


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




class diamond(fcc):
    '''Creates an instance of a diamond crystal

        Subclassed from fcc and extends that class by adding the extra atoms in
        a diamond unit cell

    '''
    def __init__(self,element,dimensionOfLattice,a):

        super().__init__(element,dimensionOfLattice,a)
        self.structure = "Diamond"

        #Hard coded diamond atoms
        extraAtoms = np.array([[0.25*a,0.25*a,0.25*a]+[0.25*a,0.75*a,0.75*a]+[0.75*a,0.75*a,0.25*a]+[0.75*a,0.25*a,0.75*a]])
        extraAtoms = extraAtoms.reshape((-1,3))

        extraAtoms = super().extendUnitCell(extraAtoms,dimensionOfLattice,a)
        self.atoms = np.concatenate((self.atoms,extraAtoms),axis = 0)




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




sc1 = sc("Si",dimensionOfLattice,LatticeConstant)

writeToxyz("sc.xyz", sc1)

bcc1 = bcc("Si",dimensionOfLattice,LatticeConstant)

writeToxyz("bcc.xyz", bcc1)

fcc1 = fcc("Si",dimensionOfLattice,LatticeConstant)

writeToxyz("fcc.xyz", fcc1)

diamond1 = diamond("Si",dimensionOfLattice,LatticeConstant)

writeToxyz("diamond.xyz", diamond1)
