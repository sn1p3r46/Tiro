import numpy as np
import scipy as sp
from scipy import linalg

from scipy.linalg.lapack import ztrsyl, dtrsyl


def sqrtm5(X,blocksize=10):
    M,Z = sp.linalg.schur(X)
    m, fb, fe, s = block_structure5(M)
    n = len(M)
    start, stop = [],[]
    start_small_blocks = []
    stop_small_blocks = []
    i, j, limit = 0, 0, 0
    while (i<n):
        start.append(i)
        limit = i + blocksize
        if (limit>n):
            limit = n
        while(i<limit):
            i+=s[j]
            j+=1
        stop.append(i)

    for i in range(0,len(start)):
        tr_sqrtm5(M[start[i]:stop[i],start[i]:stop[i]])

    for j in range(1,len(start)):
        for i in range(0,len(start)-j):
            Rii =  M[start[i]:stop[i],start[i]:stop[i]]
            Rjj =  M[start[j+i]:stop[j+i],start[j+i]:stop[j+i]]
            S = M[start[i]:stop[i],start[j+i]:stop[j+i]]
            if(j>1):
                S = S - (M[start[i]:stop[i],stop[i]:start[j+i]]).dot(
                    M[stop[i]:start[j+i],start[j+i]:stop[j+i]])

            x, scale, info = dtrsyl(Rii, Rjj, S)
            M[start[i]:stop[i],start[j+i]:stop[j+i]] = x*scale

    ZH = np.conjugate(Z).T
    X = Z.dot(M).dot(ZH)
    return X


def tr_sqrtm5(M):
    m, fb, fe, s = block_structure5(M)
    m = len(fb)
    for i in range(0,m):
        M[fb[i]:fe[i],fb[i]:fe[i]] = twobytworoot5(M[fb[i]:fe[i],fb[i]:fe[i]])
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
    #print M
    return M


def twobytworoot5(matrix):
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

def block_structure5(T):
    """
    computes the block structure of the upper quasi-triangular matrix T
    m is the number of diagonal blocks
    bb is the array containing the begin of each block
    eb is the array containing the end of each block + 1
    s is an array containing the sizes of the diagonal blocks
    """
    n = len(T)
    tol = 1e-15
    i,j = 0,0
    bb = sp.empty(n,dtype="int")
    eb = sp.empty(n,dtype="int")
    s  = sp.empty(n,dtype="int")
    while i < n-1:
        bb[j] = i
        if abs(T[i+1,i])<tol:
            i +=1
            s[j] = 1
            eb[j] = i
        else:
            i +=2
            s[j] = 2
            eb[j] = i
        j += 1
    if i == n-1:
        bb[j],eb[j] = i,i+1
        s[j] = 1
        j+=1
    bb = bb[0:j]
    eb = eb[0:j]
    s = s[0:j]
    return j, bb, eb, s
