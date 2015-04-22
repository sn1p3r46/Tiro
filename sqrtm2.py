import numpy as np
import scipy as sp
from scipy import linalg

from scipy.linalg.lapack import ztrsyl, dtrsyl


def sqrtm2(M):
    m, fb, fe, s = block_structure(M)
    n = M.shape[0]
    ris = sp.zeros((n,n))
    for i in range(0,m):
        ris[fb[i]:fe[i],fb[i]:fe[i]] = twobytworoot(M[fb[i]:fe[i],fb[i]:fe[i]])
        #print ris

    for j in range(1,m):
        for i in range(0,m-j):
            #print M[fb[i]:fe[i],fb[i+j]:fe[i+j]]
            Tnoto = sp.copy(M[fb[i]:fe[i],fb[i+j]:fe[i+j]]) #dopo togliere il copy
            #print "Tnot: "
            #print Tnoto
            for k in range(i+1,i+j):
                Tnoto -= (ris[fb[i]:fe[i],fb[k]:fe[k]]).dot(ris[fb[k]:fe[k],fb[j+i]:fe[j+i]])
                #print ris[fb[i]:fe[i],fb[k]:fe[k]]
                #print ris[fb[k]:fe[k],fb[j+i]:fe[j+i]]

            if((M[fb[i]:fe[i],fb[i+j]:fe[i+j]]).shape==(1,1)):
                #print "forma 1"
                #print M[fb[i]:fe[i],fb[i+j]:fe[i+j]]           #  Uij
                #print M[fb[i]:fe[i],fb[i]:fe[i]]               #  Uii
                #print M[fb[i+j]:fe[i+j],fb[i+j]:fe[i+j]]       #  Ujj
                ris[fb[i]:fe[i],fb[i+j]:fe[i+j]] = Tnoto/(ris[fb[i]:fe[i],fb[i]:fe[i]] + ris[fb[i+j]:fe[i+j],fb[i+j]:fe[i+j]])

            else:
                Uii = ris[fb[i]:fe[i],fb[i]:fe[i]]
                Ujj = ris[fb[i+j]:fe[i+j],fb[i+j]:fe[i+j]]
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
                x, scale, info = dtrsyl(Uii, Ujj, Tnoto)
                ris[fb[i]:fe[i],fb[i+j]:fe[i+j]] = x * scale

                """
                Tnoto = Tnoto.reshape((shapeUii*shapeUjj),1)
                ris[fb[i]:fe[i],fb[i+j]:fe[i+j]] = \
                linalg.solve(sp.kron(sp.eye(shapeUjj),Uii) +
                sp.kron(Ujj.T,sp.eye(shapeUii)),
                Tnoto).reshape(shapeUii,shapeUjj)
                """

    return ris

"""
            else:
                print "formaX"
                print M[fb[i]:fe[i],fb[i+j]:fe[i+j]]
                print M[fb[i]:fe[i],fb[i]:fe[i]]
                print M[fb[i+j]:fe[i+j],fb[i+j]:fe[i+j]]
    return ris
"""

def twobytworoot(matrix):
    if (matrix.shape[0]==2):
        a = sp.sqrt(sp.linalg.eigvals(matrix)[0]).real
        ris=sp.ndarray(shape=(2,2))
        ris[0,0] = a + (1/(4*a))*(matrix[0,0] - matrix[1,1])
        ris[1,1] = a - (1/(4*a))*(matrix[0,0] - matrix[1,1])
        ris[0,1] = (1/(2*a))*matrix[0,1]
        ris[1,0] = (1/(2*a))*matrix[1,0]
    else:
        return sp.sqrt(matrix)
    return ris

def block_structure(T):
    """
    computes the block structure of the upper quasi-triangular matrix T
    m is the number of diagonal blocks
    fb is the array containing the begin of each block
    fe is the array containing the end of each block + 1
    s is an array containing the sizes of the diagonal blocks
    """
    import scipy as sp

    n = len(T)
    tol = 1e-15

    i = 0
    v = []
    s = []

    while i < n-1:
        v.append(i)
        if abs(T[i+1,i])<tol:
            i = i + 1
            s.append(1)
        else:
            i = i + 2
            s.append(2)

    if i == n-1:
        v.append(n-1)
        s.append(1)

    m = len(v)
    fb = v
    fe = [v[i] + s[i] for i in range(m)]

    return m, fb, fe, s
