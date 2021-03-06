import numpy as np
import scipy as sp
from scipy import linalg

from scipy.linalg.lapack import ztrsyl, dtrsyl


def sqrtm4(X):
    M = sp.copy(X)
    m, fb, fe = block_structure4(M)
    n = M.shape[0]
    for i in range(0,m):
        M[fb[i]:fe[i],fb[i]:fe[i]] = twobytworoot4(M[fb[i]:fe[i],fb[i]:fe[i]])
        #print M
    for j in range(1,m):
        for i in range(0,m-j):
            JJ = i+j
            Tnoto = M[fb[i]:fe[i],fb[JJ]:fe[JJ]]
            for k in range(i+1,JJ):
                Tnoto -= (M[fb[i]:fe[i],fb[k]:fe[k]]).dot(M[fb[k]:fe[k],fb[JJ]:fe[JJ]])

            if((M[fb[i]:fe[i],fb[JJ]:fe[JJ]]).shape==(1,1)):
                M[fb[i]:fe[i],fb[JJ]:fe[JJ]] = Tnoto/(M[fb[i]:fe[i],fb[i]:fe[i]] + M[fb[JJ]:fe[JJ],fb[JJ]:fe[JJ]])

            else:
                Uii = M[fb[i]:fe[i],fb[i]:fe[i]]
                Ujj = M[fb[JJ]:fe[JJ],fb[JJ]:fe[JJ]]
                x, scale, info = dtrsyl(Uii, Ujj, Tnoto)
                M[fb[i]:fe[i],fb[JJ]:fe[JJ]] = x*scale
    return M


def twobytworoot4(matrix):
    if (matrix.shape[0]==2):
        a = sp.sqrt(sp.linalg.eigvals(matrix)[0]).real
        g = matrix[0,0] - matrix[1,1]
        matrix[0,0] = a + (1/(4*a))*(g)
        matrix[1,1] = a - (1/(4*a))*(g)
        matrix[0,1] = (1/(2*a))*matrix[0,1]
        matrix[1,0] = (1/(2*a))*matrix[1,0]
    else:
        return sp.sqrt(matrix)
    return matrix

def block_structure4(T):
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
