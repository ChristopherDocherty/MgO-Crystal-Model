#Program to produce .xyz files for Sc, Bcc, Fcc, and Diamond crystalline
#structures with the intent to display in VMD

import numpy as np


#Hard coded parameters for crystals
dimensionOfLattice = (2,2,2)
LatticeConstant = 5



class sc():
    '''

    '''

    def __init__(self,element,dimensionOfLattice,a):

        self.element = element

        self.atoms =  np.array([[0,0,0]])

        extraAtoms = self.extendUnitCell(self.atoms,dimensionOfLattice,a,(0,0,0))
        self.atoms = np.concatenate((self.atoms,extraAtoms), axis = 0)


    def extendUnitCell(self,atoms,dimensionOfLattice,a,repeatNum):
        '''

        '''

        unitVectors = [np.array([[a,0,0]]),np.array([[0,a,0]]),np.array([[0,0,a]])]

        for length, unitVector, extraTime in zip(dimensionOfLattice,unitVectors,repeatNum):
            atomCount, __  = atoms.shape

            if extraTime == 1:
                length -= 1

            for i in range(0,length):
                tempAtoms = atoms[i*atomCount:(i+1)*atomCount,:] + unitVector
                atoms = np.concatenate((atoms,tempAtoms),axis = 0)

        return atoms


class bcc(sc):

    def __init__(self,element,dimensionOfLattice,a):
        super().__init__(element,dimensionOfLattice,a)
        extraAtoms = 0.5*np.array([[a,a,a]])
        extraAtoms = super().extendUnitCell(extraAtoms,dimensionOfLattice,a,(1,1,1))
        self.atoms = np.concatenate((self.atoms,extraAtoms),axis = 0)



class fcc(sc):
    '''
    '''

    def __init__(self,element,dimensionOfLattice,a):
        super().__init__(element,dimensionOfLattice,a)
        extraAtoms = [[0,0.5*a,0.5*a],[0.5*a,0,0.5*a],[0.5*a,0.5*a,0]]
        inDim = [(0,1,1),(1,0,1),(1,1,0)]

        for atoms, repeatNum in zip(extraAtoms,inDim):

            atoms = np.array(atoms)
            atoms = atoms[np.newaxis,:]
            atoms = super().extendUnitCell(atoms,dimensionOfLattice,a,repeatNum)
            self.atoms = np.concatenate((self.atoms,atoms),axis = 0)




class diamond(fcc):
    '''

    '''
    def __init__(self,element,dimensionOfLattice,a):

        super().__init__(element,dimensionOfLattice,a)

        extraAtoms = np.array([[0.25*a,0.25*a,0.25*a]+[0.25*a,0.75*a,0.75*a]+[0.75*a,0.75*a,0.25*a]+[0.75*a,0.25*a,0.75*a]])
        extraAtoms = extraAtoms.reshape((-1,3))

        extraAtoms = super().extendUnitCell(extraAtoms,dimensionOfLattice,a,(1,1,1))
        self.atoms = np.concatenate((self.atoms,extraAtoms),axis = 0)




def writeToxyz(filename,crystal):
    '''

    '''
    crystalFile = open(filename,"wt")


    atomCountFinal = crystal.atoms.shape[0]

    #the (0,0,0) atom is always implicit
    crystalFile.write("{0}\n".format(atomCountFinal+1) )
    crystalFile.write("This is {0} \n".format(15))
    crystalFile.write("{0:s} \t 0.000000000 \t 0.000000000 \t 0.000000000 \n".format(crystal.element))
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
