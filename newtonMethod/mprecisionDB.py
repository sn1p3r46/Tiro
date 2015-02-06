import numpy as np
import scipy.linalg as spl
from mpmath import mp,workdps
import matplotlib.pyplot as plt

#mp.dps = 33
mp.dps = 128
matr = mp.matrix(
[[ '115.80200375',  '22.80402473',   '13.69453064',   '54.28049263',    '58.26628814'],
[  '22.80402473',   '86.14887381',   '53.79999432',   '42.78548627',    '37.16598947'],
[  '13.69453064',   '53.79999432',  '110.9695448' ,   '37.24270321',    '59.04033897'],
[  '54.28049263',   '42.78548627',   '37.24270321',   '95.79388469',    '27.90238338'],
[  '58.26628814',   '37.16598947',   '59.04033897',   '27.90238338',    '87.42667939']])

matr1 = np.array(
[[ 115.80200375,  22.80402473,   13.69453064,   54.28049263,  58.26628814],
[  22.80402473,   86.14887381,   53.79999432,   42.78548627, 37.16598947],
[  13.69453064,   53.79999432,  110.9695448 ,   37.24270321, 59.04033897],
[  54.28049263,   42.78548627,   37.24270321,   95.79388469, 27.90238338],
[  58.26628814,   37.16598947,   59.04033897,   27.90238338, 87.42667939]])


#Multiple Precision DB Iteration

#@workdps(128)
def MPDBiteration(M):
    Xk = M.copy()
    Yk = mp.eye(M.cols)
    for i in range(0,30):
        Xk1 = (Xk + Yk**-1)/2
        Yk1 = (Yk + Xk**-1)/2
        Xk = Xk1
        Yk = Yk1
    return Xk


D = MPDBiteration(matr)
Z = np.ndarray((matr.cols,matr.rows),dtype=np.float)

for i in range(0,M.rows):
    for j in range(0,M.cols):
        Z[i,j] = D[i,j]

def DB(M):
    Xk = np.copy(M)
    Yk = np.eye(M.shape[0])
    res = []
    for i in range(0,30):
        Xk1 = (Xk + np.linalg.inv(Yk))/2
        Yk1 = (Yk + np.linalg.inv(Xk))/2
        Xk = Xk1
        Yk = Yk1
        res.append(spl.norm(Xk-Z)/spl.norm(Z))

    return Xk,res

def newton(M):
    A = np.copy(M)
    X = np.copy(M)
    res = []
    for i in range(0,30):
        X = 0.5*(X + np.linalg.inv(X).dot(A))
        res.append(spl.norm(X-Z)/spl.norm(Z))

    return X,res

def productDB(M):
    Xk = np.copy(M)
    Mk = np.copy(M)
    I = np.eye(M.shape[0])
    res = []
    for i in range(0,30):
        Xk = (0.5)*(Xk).dot(I + np.linalg.inv(Mk))
        Mk = (I + ((Mk + np.linalg.inv(Mk))/2))/2
        res.append(spl.norm(Xk-Z)/spl.norm(Z))

    return Xk,res

def CRiteration(M):
    I = np.eye(M.shape[0])
    Yk = I - M
    Zk = 2*(I + M)
    res = []
    for i in range(0,30):
        Yk = ((-Yk.dot(np.linalg.inv(Zk)).dot(Yk)))
        Zk = Zk + 2 * Yk
        res.append(spl.norm(Zk/4-Z)/spl.norm(Z))

    return Zk,res

def INiteration(M):
    I = np.eye(M.shape[0])
    Xk = np.copy(M)
    Ek = (I - M)/2
    res = []
    for i in range(0,30):
        Xk = Xk + Ek
        Ek = ((Ek.dot(np.linalg.inv(Xk)).dot(Ek)))/(-2)
        res.append(spl.norm(Xk-Z)/spl.norm(Z))

    return Xk,res

a = np.arange(0,30)

nw = [0]*2
db = [0]*2
dbpr = [0]*2
iN = [0]*2

nw[0],nw[1] = newton(matr1)
db[0],db[1] = DB(matr1)
dbpr[0],dbpr[1] = productDB(matr1)
iN[0],iN[1] = INiteration(matr1)

plt.plot(a,iN[1])
plt.plot(a,dbpr[1])
plt.plot(a,db[1])
plt.plot(a,nw[1])
plt.semilogy(t, np.exp(-t/5.0))
plt.grid(True)
plt.show()


#[(x,y) for x in range(0,3) for y in range(0,3)]
