import numpy as np
import scipy as sp
import scipy.linalg
import random


def complexSchur(A,array=None):
    if array is None:
        array = 1
    Res = sp.diag(A)*1.0
    Res = sp.diag(Res)
    for i in range(0,len(Res)):
        Res[i] = sp.sqrt(Res[i])
    Res = Res*array
    for j in range(1,A.shape[0]):
        for i in range(0,A.shape[0]-j):
            s = 0
            tij = A[i,j+i]
            den = Res[i,i] + Res[j+i,j+i]
            for k in range(i+1,j+i):
                s += (Res[i,k]*Res[k,j+i])
            Res[i,j+i] = (tij-s)/(den)
    return Res


if __name__ == "__main__":
    #A = sp.array([[  29.15568676-0.j,   41.00000000+0.j],[   0.00000000+0.j,  125.15568676+0.j]])

    A = sp.array([[8.-2j, 93., 28., 36.],
                  [0.,  5., 14., 27.],
                  [0.,  0.,  4., 86.],
                  [0.,  0.,  0.,  2.]])

    matrice = complexSchur(A)
    print matrice.dot(matrice)
    print sp.linalg.norm(matrice.dot(matrice)-A)
    print (matrice.dot(matrice)-A).max()
