import scipy as sp
import numpy as np

def newton(M):
    A = np.copy(M)
    X = np.copy(M)
    for i in range(0,30):
        X = 0.5*(X + np.linalg.inv(X).dot(A))

    return X

def DB(M):
    Xk = np.copy(M)
    Yk = np.eye(M.shape[0])
    for i in range(0,30):
        Xk1 = (Xk + np.linalg.inv(Yk))/2
        Yk1 = (Yk + np.linalg.inv(Xk))/2
        Xk = Xk1
        Yk = Yk1
    return Xk

def productDB(M):
    Xk = np.copy(M)
    Mk = np.copy(M)
    I = np.eye(M.shape[0])
    for i in range(0,30):
        Xk = (0.5)*(Xk).dot(I + np.linalg.inv(Mk))
        Mk = (I + ((Mk + np.linalg.inv(Mk))/2))/2
    return Xk

def CRiteration(M):
    I = np.eye(M.shape[0])
    Yk = I - M
    Zk = 2*(I + M)
    for i in range(0,30):
        Yk = ((-Yk.dot(np.linalg.inv(Zk)).dot(Yk)))
        Zk = Zk + 2 * Yk

    return Zk

def INiteration(M):
    I = np.eye(M.shape[0])
    Xk = np.copy(M)
    Ek = (I - M)/2
    for i in range(0,30):
        Xk = Xk + Ek
        Ek = ((Ek.dot(np.linalg.inv(Xk)).dot(Ek)))/(-2)
    return Xk



"""

[[87, 87, 15, 14, 70],
[22, 47, 89, 53, 25],
[62, 16, 82, 97,  6],
[74, 73, 72,  7,  9],
[64, 58, 30, 82, 40]]

"""
"""

[[31, 95, 26, 78],
[25, 29, 86, 29],
[97, 48, 46, 49],
[ 7, 27, 45, 45]]

"""
