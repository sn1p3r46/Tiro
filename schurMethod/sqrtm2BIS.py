import numpy as np
import scipy as sp
from scipy import linalg
from scipy.linalg.lapack import ztrsyl, dtrsyl
from os.path import expanduser
import os
import scipy as sp
import scipy.linalg
import time



execfile("myHand.py")
execfile("aggiustala.py")

def realSchurNonBlocked(Matrix):
    M,T = sp.linalg.schur(Matrix)
    m, fb, fe, s = block_structure(M) #computes the block Structure
    n = M.shape[0]
    ris = sp.zeros((n,n))
    for i in range(0,m):
        ris[fb[i]:fe[i],fb[i]:fe[i]] = tbtroot(M[fb[i]:fe[i],fb[i]:fe[i]])

    for j in range(1,m):
        for i in range(0,m-j):
            Tnoto = sp.copy(M[fb[i]:fe[i],fb[i+j]:fe[i+j]]) #dopo togliere il copy

            for k in range(i+1,i+j):
                Tnoto -= (ris[fb[i]:fe[i],fb[k]:fe[k]]).dot(ris[fb[k]:fe[k],fb[j+i]:fe[j+i]])

            if((M[fb[i]:fe[i],fb[i+j]:fe[i+j]]).shape==(1,1)):
                ris[fb[i]:fe[i],fb[i+j]:fe[i+j]] = Tnoto/(ris[fb[i]:fe[i],fb[i]:fe[i]] + ris[fb[i+j]:fe[i+j],fb[i+j]:fe[i+j]])

            else:
                Uii = ris[fb[i]:fe[i],fb[i]:fe[i]]
                Ujj = ris[fb[i+j]:fe[i+j],fb[i+j]:fe[i+j]]
                shapeUii = Uii.shape[0]
                shapeUjj = Ujj.shape[0]
                x, scale, info = dtrsyl(Uii, Ujj, Tnoto) #call fortran and solve the slvester equation
                ris[fb[i]:fe[i],fb[i+j]:fe[i+j]] = x * scale
    return ris


def csr2(complex_n):
    t = sp.sqrt((sp.absolute(complex_n.real)+(sp.sqrt(sp.square(complex_n.real)+sp.square(complex_n.imag))))/2)
    if(complex_n.real >= 0):
        return complex(t,(complex_n.imag/(2*t)))
    else:
        return complex((complex_n.imag/(2*t)),t)

def gev(matrix):
    x = (matrix[0,0]+matrix[1,1])/2
    y = sp.sqrt(-sp.square(matrix[0,0]-matrix[1,1])-4*matrix[1,0]*matrix[0,1])/2
    return complex(x,y)

def tbtroot(matrix):
    if (matrix.shape[0]==2):
        a = gev(matrix)
        a = csr2(a).real
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


def timeTestRealSchur(M):
    I = time.clock()
    Res = realSchurNonBlocked(M)
    F = time.clock()
    DELTA = F - I
    ERR = sp.linalg.norm(Res.dot(Res) - M)/sp.linalg.norm(M)

    return DELTA,ERR

B = scipy.random.randint(-100,100,(10,10))
print "\nCreo matrice B"
print "\nOra calcolo schur di B!"
T,Z = sp.linalg.schur(B)
print "\nDecomposizione di schur terminata con sucesso !!!\n "
aggiustala(T)

DD = realSchurNonBlocked(T)

print sp.linalg.norm(DD.dot(DD)-T)/sp.linalg.norm(T)
