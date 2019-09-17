#Program to produce .xyz files for Sc, Bcc, Fcc, and Diamond crystalline
#structures with the intent to display in VMD

import numpy as np


#Hard coded parameters for crystals
dimensionOfLattice = (1,1,1)
LatticeConstant = 5



class sc():
    '''

    '''
    extraAtoms = np.array([])


    def __init__(self,structure,element,dimensionOfLattice,a):

        self.structure = structure
        self.element = element

        self.atoms =  np.array([[0,0,0]])
        self.atoms = self.extendLatticePoint(self.atoms,dimensionOfLattice,a,0)


    def extendLatticePoint(self,atoms,dimensionOfLattice,a,nonVertex):
        '''
        '''

        def unitCellGenerator(direction):
            '''
            '''

            length = dimensionOfLattice[direction] - nonVertex


            #atomCount is the number in a unit cell
            atomCount, __  = atoms.shape

            #Need to RENAME
            if direction == 0:
                vector = np.array([[a,0,0]])
            elif direction == 1:
                vector = np.array([[0,a,0]])
            elif direction == 2:
                vector = np.array([[0,0,a]])

            i = 0
            #SO MESSY
            while i < length:

                #Dimensions of atomsInBlock x (3*#atoms)
                newUnitCells = atoms[i*atomCount:(i+1)*atomCount,:] + vector
                i += 1
                #probably unecessary to keep converting - just one at the end

                yield newUnitCells



        #size among dimension is stored in a tuple so 0 accesses xlenght etc.


        xExtensionGen = unitCellGenerator(0)
        for k in xExtensionGen:
            atoms = np.concatenate((atoms,k), axis=0)
        yExtensionGen = unitCellGenerator(1)
        for k in yExtensionGen:
            atoms = np.concatenate((atoms,k), axis=0)
        zExtensionGen = unitCellGenerator(2)
        for k in zExtensionGen:
            atoms = np.concatenate((atoms,k), axis=0)

        return atoms




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


class bcc(sc):

    def __init__(self,structure,element,dimensionOfLattice,a):
        super(bcc, self).__init__(structure,element,dimensionOfLattice,a)
        extraAtoms = 0.5*np.array([[a,a,a]])
        extraAtoms = self.extendLatticePoint(extraAtoms,dimensionOfLattice,a,1)
        self.atoms = np.concatenate((self.atoms,extraAtoms),axis = 0)


class fcc(sc):

    def __init__(self,structure,element,dimensionOfLattice,a):
        super(fcc, self).__init__(structure,element,dimensionOfLattice,a)
        a2 = a/2
        extraAtoms =np.array([[a2,a2,0]+[a2,0,a2]+[0,a2,a2]])
        extraAtoms = extraAtoms.reshape((-1,3))
        self.extendLatticePoint(extraAtoms,dimensionOfLattice,a,0)
        self.atoms = np.concatenate((self.atoms,extraAtoms),axis = 0)



sc = sc("sc","Si",dimensionOfLattice,LatticeConstant)

writeToxyz("sc.xyz",sc)

bcc = bcc("bcc","Si",dimensionOfLattice,LatticeConstant)

writeToxyz("bcc.xyz",bcc)

fcc = fcc("fcc","Si",dimensionOfLattice,LatticeConstant)

writeToxyz("fcc.xyz",fcc)

#diamond = crystal("diamond","Si",dimensionOfLattice,LatticeConstant)

#writeToxyz("diamond.xyz",diamond)
