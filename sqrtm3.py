import numpy as np
import scipy as sp
from scipy import linalg

from scipy.linalg.lapack import ztrsyl, dtrsyl


def sqrtm3(X):
    M = sp.copy(X)
    m, fb, fe = block_structure(M)
    n = M.shape[0]
    for i in range(0,m):
        M[fb[i]:fe[i],fb[i]:fe[i]] = twobytworoot(M[fb[i]:fe[i],fb[i]:fe[i]])
        #print M

    for j in range(1,m):
        for i in range(0,m-j):
            #print M[fb[i]:fe[i],fb[JJ]:fe[JJ]]
            JJ = i+j
            Tnoto = M[fb[i]:fe[i],fb[JJ]:fe[JJ]] #dopo togliere il copy
            #print "Tnot: "
            #print Tnoto
            for k in range(i+1,JJ):
                Tnoto -= (M[fb[i]:fe[i],fb[k]:fe[k]]).dot(M[fb[k]:fe[k],fb[JJ]:fe[JJ]])
                #print M[fb[i]:fe[i],fb[k]:fe[k]]
                #print M[fb[k]:fe[k],fb[JJ]:fe[JJ]]

            if((M[fb[i]:fe[i],fb[JJ]:fe[JJ]]).shape==(1,1)):
                #print "forma 1"
                #print M[fb[i]:fe[i],fb[JJ]:fe[JJ]]           #  Uij
                #print M[fb[i]:fe[i],fb[i]:fe[i]]               #  Uii
                #print M[fb[JJ]:fe[JJ],fb[JJ]:fe[JJ]]       #  Ujj
                M[fb[i]:fe[i],fb[JJ]:fe[JJ]] = Tnoto/(M[fb[i]:fe[i],fb[i]:fe[i]] + M[fb[JJ]:fe[JJ],fb[JJ]:fe[JJ]])

            else:
                Uii = M[fb[i]:fe[i],fb[i]:fe[i]]
                Ujj = M[fb[JJ]:fe[JJ],fb[JJ]:fe[JJ]]
                shapeUii = Uii.shape[0]
                shapeUjj = Ujj.shape[0]
                """
                print "------------"
                print Tnoto
                print Tnoto.shape
                print sp.kron(sp.eye(shapeUjj),Uii)
                print sp.kron(Ujj.T,sp.eye(shapeUii))
                print Tnoto
                """
                #M[fb[i]:fe[i],fb[JJ]:fe[JJ]] = sp.linalg.solve_sylvester(Uii, Ujj, Tnoto)

                """
                x, scale, info = dtrsyl(Uii, Ujj, Tnoto

                if (scale==1.0):
                     = x

                else:
                    M[fb[i]:fe[i],fb[JJ]:fe[JJ]] = x*scale
                    print "scale!=0"
                """
                Tnoto = Tnoto.reshape((shapeUii*shapeUjj),1,order="F")
                M[fb[i]:fe[i],fb[JJ]:fe[JJ]] = \
                linalg.solve(sp.kron(sp.eye(shapeUjj),Uii) +
                sp.kron(Ujj.T,sp.eye(shapeUii)),
                Tnoto).reshape(shapeUii,shapeUjj,order="F")


    return M

"""
            else:
                print "formaX"
                print M[fb[i]:fe[i],fb[JJ]:fe[JJ]]
                print M[fb[i]:fe[i],fb[i]:fe[i]]
                print M[fb[JJ]:fe[JJ],fb[JJ]:fe[JJ]]
    return M
"""

def twobytworoot(matrix):
    if (matrix.shape[0]==2):
        a = sp.sqrt(sp.linalg.eigvals(matrix)[0]).real
        M=sp.ndarray(shape=(2,2))
        M[0,0] = a + (1/(4*a))*(matrix[0,0] - matrix[1,1])
        M[1,1] = a - (1/(4*a))*(matrix[0,0] - matrix[1,1])
        M[0,1] = (1/(2*a))*matrix[0,1]
        M[1,0] = (1/(2*a))*matrix[1,0]
    else:
        return sp.sqrt(matrix)
    return M

def block_structure(T):
    """
    computes the block structure of the upper quasi-triangular matrix T
    m is the number of diagonal blocks
    fb is the array containing the begin of each block
    fe is the array containing the end of each block + 1
    s is an array containing the sizes of the diagonal blocks
    """
    n = len(T)
    tol = 1e-15
    i,j = 0,0
    bb = sp.empty(n,dtype="int")
    eb = sp.empty(n,dtype="int")
    while i < n-1:
        bb[j] = i
        if abs(T[i+1,i])<tol:
            i +=1
            eb[j] = i
        else:
            i +=2
            eb[j] = i
        j += 1
    if i == n-1:
        bb[j],eb[j] = i,i+1
        j+=1
    return j, bb, eb
