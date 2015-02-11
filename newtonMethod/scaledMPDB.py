import numpy as np
import scipy.linalg as spl
from mpmath import mp,workdps
import matplotlib.pyplot as plt
from mprecisionDB import DB
from mprecisionDB import MPDBiteration,matr,matr1

class result:
    def __init__(self,iter=0,ris=None,res=[]):
        self.iter = iter
        self.ris = ris
        self.res = res

def MPDBsc(M):
    Xk = M.copy()
    Yk = mp.eye(M.cols)
    dA = mp.det(M)**(mp.mpf('1.0')/2)
    for i in range(1,30):
        uk = mp.fabs(mp.det(Xk)/dA)**(-1.0/i)
        Xk1 = (uk*Xk + (uk**(-1))*(Yk**(-1)))/2
        Yk1 = (uk*Yk + (uk**(-1))*(Xk**(-1)))/2
        Xk = Xk1
        Yk = Yk1

    return Xk


def NWsc(M, object=False):
    W = np.array(MPDBsc((mp.matrix(M))).tolist(),dtype=np.float64)
    X = np.copy(M)
    dM = spl.det(M)**(1.0/2)
    res = result()
    for i in range(1,31):
        uk = np.abs(spl.det(X)/dM)**(-1.0/i)
        X = (0.5)*(uk*X + (uk**(-1))*spl.inv(X).dot(M))
        res.res.append(spl.norm(X-W)/spl.norm(W))
    res.iter = i
    res.ris = X

    return (res if object else X)


def DBsc(M, object=False):
    W = np.array(MPDBsc((mp.matrix(M))).tolist(),dtype=np.float64)
    res = result()
    Xk = np.copy(M)
    Yk = np.eye(M.shape[0])
    dA = spl.det(M)**(1.0/2)
    for i in range(1,31):
        uk = np.abs(spl.det(Xk)/dA)**(-1.0/i)
        Xk1 = (Xk*uk + (uk**-1)*np.linalg.inv(Yk))/2
        Yk1 = (Yk*uk + (uk**-1)*np.linalg.inv(Xk))/2
        Xk = Xk1
        Yk = Yk1
        res.res.append(spl.norm(Xk-W)/spl.norm(W))
    res.iter = i
    res.ris = Xk

    return (res if object else Xk)

def productDBsc(A, object=False):
    W = np.array(MPDBsc((mp.matrix(A))).tolist(),dtype=np.float64)
    res = result()
    M = np.copy(A)
    Xk = np.copy(A)
    Yk = np.eye(A.shape[0])
    I = np.eye(M.shape[0])
    for i in range(1,31):
        uk = np.abs(spl.det(M))**(-1.0/(2*i))
        Xk = (0.5)*(uk)*(Xk).dot(I + ((uk**-2)*(np.linalg.inv(M))))
        Yk = (0.5)*(uk)*(Yk).dot(I + ((uk**-2)*(np.linalg.inv(M))))
        M  = (0.5)*( I + ((uk**2)*M + ((uk**-2)*np.linalg.inv(M)))/2.0)
        res.res.append(spl.norm(Xk-W)/spl.norm(W))
    res.iter = i
    res.ris = Xk

    return (res if object else Xk)

"""
def INsc1(A):
    X = np.copy(A)
    E = 0.5*(np.eye(A.shape[0])-A)
    dA = spl.det(A)**(0.5)
    for i in range (1,31):
        print i
        detX =  spl.det(X)
        uk = np.abs(detX/dA)**(-1.0/i)
        Etk = (uk**(-1))*(E + (0.5*X)) - (0.5)*uk*X
        X = uk*X + Etk
        E = (-0.5)*(E.dot(spl.inv(X))).dot(E)
    return X
"""
