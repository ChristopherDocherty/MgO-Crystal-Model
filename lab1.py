#Program to produce .xyz files for Sc, Bcc, Fcc, and Diamond crystalline
#structures with the intent to display in VMD

import numpy as np


#Hard coded parameters for crystals
dimensionOfLattice = (1,1,1)
LatticeConstant = 5


class crystal():
    '''

    '''

    a = LatticeConstant
    def __init__(self,structure,element,dimensionOfLattice,LatticeConstant):

        self.structure = structure
        self.element = element
        a = LatticeConstant
        #Each column contains a unit cell
        #The base atom for the unti cell is the last of the previous column
        #except for the first where it is always assumed to be (0,0,0)
        self.atoms =  np.array([[a,0,0]+[0,a,0]+[0,0,a]+[a,a,0]+[a,0,a]+[a,a,a]+[0,a,a]])
        self.atoms = self.atoms.reshape((-1,3))
        self.addExtraAtoms(structure,a)

        self.extendUnitCell(dimensionOfLattice,a)



    def addExtraAtoms(self,structure,a):
        '''
        '''
        if structure == "sc":
            return
        elif structure == "bcc":
            self.atoms =np.concatenate((self.atoms,0.5*np.array([[a,a,a]])),axis=0)
        elif structure == "fcc" or structure == "diamond":
            a2 = 0.5*a

            extraAtoms =np.array([[a2,a2,0]+[a2,0,a2]+[0,a2,a2]+[a,a2,a2]+[a2,a2,a]+[a2,a,a2]])
            extraAtoms = extraAtoms.reshape((-1,3))
            self.atoms = np.concatenate((self.atoms,extraAtoms),axis=0)

            if structure == "diamond":
                extraAtoms = np.array([[0.25*a,0.25*a,0.25*a]+[0.25*a,0.25*a,0.75*a]+[0.5*a,0.75*a,0.25*a]+[0.75*a,0.75*a,0.75*a]])
                extraAtoms = extraAtoms.reshape((-1,3))
                self.atoms = np.concatenate((self.atoms,extraAtoms),axis=0)

        return


    def extendUnitCell(self,dimensionOfLattice,LatticeConstant):
        '''
        '''

        def unitCellGenerator(direction):
            '''
            '''

            length = dimensionOfLattice[direction]


            #atomCount is the number in a unit cell
            atomCount, __  = self.atoms.shape

            #Need to RENAME
            if direction == 0:
                vector = np.array([[LatticeConstant,0,0]])
            elif direction == 1:
                vector = np.array([[0,LatticeConstant,0]])
            elif direction == 2:
                vector = np.array([[0,0,LatticeConstant]])

            i = 0
            #SO MESSY
            while i < length-1:

                #Dimensions of atomsInBlock x (3*#atoms)
                newUnitCells = self.atoms[i*atomCount:(i+1)*atomCount,:] + vector
                i += 1
                #probably unecessary to keep converting - just one at the end

                yield newUnitCells



        #size among dimension is stored in a tuple so 0 accesses xlenght etc.


        xExtensionGen = unitCellGenerator(0)
        for k in xExtensionGen:
            self.atoms = np.concatenate((self.atoms,k), axis=0)


        yExtensionGen = unitCellGenerator(1)
        for k in yExtensionGen:
            self.atoms = np.concatenate((self.atoms,k), axis=0)

        zExtensionGen = unitCellGenerator(2)
        for k in zExtensionGen:
            self.atoms = np.concatenate((self.atoms,k), axis=0)

        return


#Need to add file writing function
def writeToxyz(filename,crystal):
    '''

    '''
    crystalFile = open(filename,"wt")


    atomCountFinal = crystal.atoms.shape[0]

    #the (0,0,0) atom is always implicit
    crystalFile.write("{0}\n".format(atomCountFinal+1) )
    crystalFile.write("This is {0} \n".format(sc.structure))
    crystalFile.write("{0:s} \t 0.000000000 \t 0.000000000 \t 0.000000000 \n".format(crystal.element))
    for i in range(0,atomCountFinal):
        crystalFile.write("{0:s} \t {1:0.9f} \t {2:0.9f} \t {3:0.9f} \n".format(crystal.element,crystal.atoms[i,0],crystal.atoms[i,1],crystal.atoms[i,2]))




    return


sc = crystal("sc","Si",dimensionOfLattice,LatticeConstant)

writeToxyz("sc.xyz",sc)

bcc = crystal("bcc","Si",dimensionOfLattice,LatticeConstant)

writeToxyz("bcc.xyz",bcc)

fcc = crystal("fcc","Si",dimensionOfLattice,LatticeConstant)

writeToxyz("fcc.xyz",fcc)

diamond = crystal("diamond","Si",dimensionOfLattice,LatticeConstant)

writeToxyz("diamond.xyz",diamond)
