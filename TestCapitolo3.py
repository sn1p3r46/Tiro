import scipy as sp
import numpy as np
import scipy.linalg as spl

import matplotlib
import matplotlib.pyplot as plt

moler = np.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/moler.txt")
frank = np.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/frank.txt")
vand  = np.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/vand.txt")

matr1 = np.array(
[[ 115.80200375,  22.80402473,   13.69453064,   54.28049263,  58.26628814],
[  22.80402473,   86.14887381,   53.79999432,   42.78548627, 37.16598947],
[  13.69453064,   53.79999432,  110.9695448 ,   37.24270321, 59.04033897],
[  54.28049263,   42.78548627,   37.24270321,   95.79388469, 27.90238338],
[  58.26628814,   37.16598947,   59.04033897,   27.90238338, 87.42667939]])



class result:
    def __init__(self,iter=0,ris=None, res=[]):
        self.iter = 0
        self.ris = 0
        self.res = []

def NWsc(M, object=False):
    X = np.copy(M)
    dM = spl.det(M)**(1.0/2)
    ras = result()
    for i in range(1,31):
        print i
        uk = np.abs(spl.det(X)/dM)**(-1.0/i)
        X = (0.5)*(uk*X + (uk**(-1))*spl.inv(X).dot(M))
        print spl.norm(X.dot(X)-M)/spl.norm(M)
        ras.res.append((spl.norm(X.dot(X)-M)/spl.norm(M)))
    ras.iter = i
    ras.ris = X
    return (ras if object else X)

def DBsc(M, object=False):
    ras = result()
    Xk = np.copy(M)
    Yk = np.eye(M.shape[0])
    dA = spl.det(M)**(1.0/2)
    for i in range(1,31):
        uk = np.abs(spl.det(Xk)/dA)**(-1.0/i)
        Xk1 = (Xk*uk + (uk**-1)*np.linalg.inv(Yk))/2
        Yk1 = (Yk*uk + (uk**-1)*np.linalg.inv(Xk))/2
        Xk = Xk1
        Yk = Yk1
        ras.res.append((spl.norm(Xk.dot(Xk)-M)/spl.norm(M)))
    ras.iter = i
    ras.ris = Xk

    return (ras if object else Xk)

def productDBsc(A, object=False):
    ras = result()
    M = np.copy(A)
    Xk = np.copy(A)
    Yk = np.eye(A.shape[0])
    I = np.eye(M.shape[0])
    for i in range(1,31):
        uk = np.abs(spl.det(M))**(-1.0/(2*i))
        Xk = (0.5)*(uk)*(Xk).dot(I + ((uk**-2)*(np.linalg.inv(M))))
        Yk = (0.5)*(uk)*(Yk).dot(I + ((uk**-2)*(np.linalg.inv(M))))
        M  = (0.5)*( I + ((uk**2)*M + ((uk**-2)*np.linalg.inv(M)))/2.0)
        ras.res.append((spl.norm(Xk.dot(Xk)-A)/spl.norm(A)))
    ras.iter = i
    ras.ris = Xk

    return (ras if object else Xk)

def INsc1(A, object=False):
    X = np.copy(A)
    E = 0.5*(np.eye(A.shape[0])-A)
    dA = spl.det(A)**(0.5)
    ras = result()
    for i in range (1,31):
        print i
        detX =  spl.det(X)
        uk = np.abs(detX/dA)**(-1.0/i)
        Etk = (uk**(-1))*(E + (0.5*X)) - (0.5)*uk*X
        X = uk*X + Etk
        E = (-0.5)*(Etk.dot(spl.inv(X))).dot(Etk)
        ras.res.append((spl.norm(X.dot(X)-A)/spl.norm(A)))
        print ras.res
    ras.iter = i
    ras.ris = X
    return (ras if object else X)

def plottamelo():
    a = np.arange(0,30)
    nw = NWsc(matr1,True)
    db = DBsc(matr1,True)
    dbPR = productDBsc(matr1,True)
    INsc = INsc1(matr1,True)
    plt.yscale('log')
    plt.plot(range(0,nw.iter),nw.res, color='blue', lw=2, label="Newton Sclaled")
    plt.plot(range(0,db.iter),db.res, color='red', lw=2, label="DB Sclaled")
    plt.plot(range(0,dbPR.iter),dbPR.res, color='green', lw=2, label="DB Product Sclaled")
    plt.plot(range(0,INsc.iter),INsc.res, color='orange', lw=2, label="IN Sclaled")
    plt.legend(loc='upper right')
    plt.show()



def newton(M,object=False):
    ras = result()
    A = np.copy(M)
    X = np.copy(M)
    for i in range(0,30):
        X = 0.5*(X + np.linalg.inv(X).dot(A))
        ras.res.append((spl.norm(X.dot(X)-M)/spl.norm(M)))
    ras.iter = i
    ras.ris = X
    return (ras if object else X)

def DB(M,object=False):
    ras = result()
    Xk = np.copy(M)
    Yk = np.eye(M.shape[0])
    for i in range(0,30):
        Xk1 = (Xk + np.linalg.inv(Yk))/2
        Yk1 = (Yk + np.linalg.inv(Xk))/2
        Xk = Xk1
        Yk = Yk1
        ras.res.append((spl.norm(Xk.dot(Xk)-M)/spl.norm(M)))
    ras.iter = i
    ras.ris = Xk
    return (ras if object else Xk)

def productDB(M,object=False):
    ras = result()
    Xk = np.copy(M)
    Mk = np.copy(M)
    I = np.eye(M.shape[0])
    for i in range(0,30):
        Xk = (0.5)*(Xk).dot(I + np.linalg.inv(Mk))
        Mk = (I + ((Mk + np.linalg.inv(Mk))/2))/2
        ras.res.append((spl.norm(Xk.dot(Xk)-M)/spl.norm(M)))
    ras.iter = i
    ras.ris = Xk
    return (ras if object else Xk)

def CRiteration(M,object=False):
    ras = result()
    I = np.eye(M.shape[0])
    Yk = I - M
    Zk = 2*(I + M)
    for i in range(0,30):
        Yk = ((-Yk.dot(np.linalg.inv(Zk)).dot(Yk)))
        Zk = Zk + 2 * Yk
        ras.res.append((spl.norm((Zk/4.0).dot((Zk/4.0))-M)/spl.norm(M)))
    ras.iter = i
    ras.ris = Zk/4.0
    return (ras if object else Zk)

def INiteration(M,object=False):
    ras = result()
    I = np.eye(M.shape[0])
    Xk = np.copy(M)
    Ek = (I - M)/2
    for i in range(0,30):
        Xk = Xk + Ek
        Ek = ((Ek.dot(np.linalg.inv(Xk)).dot(Ek)))/(-2)
        ras.res.append((spl.norm(Xk.dot(Xk)-M)/spl.norm(M)))
    ras.iter = i
    ras.ris = Xk
    return (ras if object else Xk)


def plottamelo2():
    a = np.arange(0,30)
    nw = newton(matr1,True)
    db = DB(matr1,True)
    dbPR = productDB(matr1,True)
    IN = INiteration(matr1,True)
    CR = CRiteration(matr1,True)
    plt.yscale('log')
    print len(range(0,nw.iter))
    print len(nw.res)
    plt.plot(range(0,nw.iter+1),nw.res, color='blue', lw=2, label="Newton")
    plt.plot(range(0,db.iter+1),db.res, color='red', lw=2, label="DB")
    plt.plot(range(0,dbPR.iter+1),dbPR.res, color='green', lw=2, label="DB Product")
    plt.plot(range(0,IN.iter+1),IN.res, color='orange', lw=2, label="IN")
    plt.plot(range(0,CR.iter+1),CR.res, color='brown', lw=2, label="CR")

    plt.legend(loc='upper right')
    plt.show()


def plottamelo3():
    a = np.arange(0,30)
    nw = NWsc(matr1,True)
    db = DBsc(matr1,True)
    dbPR = productDBsc(matr1,True)
    INsc = INsc1(matr1,True)
    plt.yscale('log')
    plt.plot(range(0,nw.iter),nw.res, color='blue', lw=2, label="Newton Sclaled")
    #plt.plot(range(0,db.iter),db.res, color='red', lw=2, label="DB Sclaled")
    #plt.plot(range(0,dbPR.iter),dbPR.res, color='green', lw=2, label="DB Product Sclaled")
    #plt.plot(range(0,INsc.iter),INsc.res, color='orange', lw=2, label="IN Sclaled")
    nw = newton(matr1,True)
    db = DB(matr1,True)
    dbPR = productDB(matr1,True)
    IN = INiteration(matr1,True)
    CR = CRiteration(matr1,True)
    plt.yscale('log')
    print len(range(0,nw.iter))
    print len(nw.res)
    plt.plot(range(0,nw.iter+1),nw.res, color='blue', lw=2, label="Newton")
    #plt.plot(range(0,db.iter+1),db.res, color='red', lw=2, label="DB")
    #plt.plot(range(0,dbPR.iter+1),dbPR.res, color='green', lw=2, label="DB Product")
    #plt.plot(range(0,IN.iter+1),IN.res, color='orange', lw=2, label="IN")
    #plt.plot(range(0,CR.iter+1),CR.res, color='brown', lw=2, label="CR")

    plt.legend(loc='upper right')
    plt.show()
