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



def steepestDescent(g, f, pos_vec, a=0.1, h = 10**(-8)):



    g0 = g(pos_vec)

    delg = f(pos_vec)

    print(delg)
    dg = math.sqrt(dotProduct(delg,delg))

    b = a/dg

    count = 0
    while dg > h:



        pos_vec -= b * delg

        delg = f(pos_vec)
        dg =  math.sqrt(dotProduct(delg,delg))
        b = a/dg



        g1 = g(pos_vec)
        if g1 > g0:
            a /= 2
        else:
            g0 = g1

        count += 1
    print(count)
    return pos_vec


x = np.array([2.0,3.0,4.0])



result =steepestDescent(LJ_potential,LJ_derivative,x)

print(2**(1/6) * sigma)
print(result)
print(math.sqrt(dotProduct(result,result)))
