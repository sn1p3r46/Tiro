import numpy as np
import scipy as sp
from scipy import linalg
from scipy.linalg.lapack import ztrsyl, dtrsyl
from os.path import expanduser
import os
import scipy as sp
import scipy.linalg
import time
from os.path import expanduser



workspace = "/Tiro/IOFiles/"
folder = "/TestCapitolo2/"
home = expanduser("~")

execfile("schurMethod/sqrtm2BIS.py")
execfile("schurMethod/complexSchur.py")


MATRICE = sp.loadtxt("IOFiles/BigBlockMatrix5000.out")

def timeTestCap2(MATRICE,dim):
    Dim = 10
    TimeReal=[]
    TimeComplex =[]
    ErrReal = []
    ErrComplex =[]
    DIMS = []
    for i in range(0,dim):
        print Di
        DIMS.append(Dim)
        M = MATRICE[0:Dim,0:Dim]

        I = time.clock()
        Res = realSchurNonBlocked(M)
        F = time.clock()
        TimeReal.append(F - I)
        ErrReal.append(sp.linalg.norm(Res.dot(Res) - M)/sp.linalg.norm(M))

        I = time.clock()
        Res = complexSchur(M)
        F = time.clock()
        TimeComplex.append(F - I)
        ErrComplex.append(sp.linalg.norm(Res.dot(Res) - M)/sp.linalg.norm(M))
        Dim = Dim*2

    print ErrReal
    print ErrComplex
    print TimeReal
    print TimeComplex

    np.savetxt(home + workspace + folder + "TimeReal"+str(dim)+".tt", TimeReal)
    np.savetxt(home + workspace + folder + "TimeComplex"+str(dim)+".tt", TimeComplex)
    np.savetxt(home + workspace + folder + "ErrReal"+str(dim)+".tt", ErrReal)
    np.savetxt(home + workspace + folder + "ErrComplex"+str(dim)+".tt", ErrComplex)
    np.savetxt(home + workspace + folder + "DIMS"+str(dim)+".tt", DIMS)




def rarra(cicli):
    Dim = 10
    for i in  range(0,cicli):
        print Dim
        Dim = Dim*2
