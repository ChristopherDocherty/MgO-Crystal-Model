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





def steepestDescent(fprime, pos_vec, h = 10**(-7)):


    dg = 1
    count = 0

    for i in range(0,13):
        g = -fprime(pos_vec)
        dg = math.sqrt(dotProduct(g,g))

        a = aGet(pos_vec,fprime,g)
        #print("pos_vec = {0}   g = {1}   a = {2} hmm {3}".format(pos_vec,g,a,g*a))
        pos_vec += a * g
        count+= 1

    #print(count)
    return pos_vec


x = np.array([1.0,1.0,1.0])



result =steepestDescent(LJ_derivative,x)

print(2**(1/6) * sigma)
print(result)
print(math.sqrt(dotProduct(result,result)))
