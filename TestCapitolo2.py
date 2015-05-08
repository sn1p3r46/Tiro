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
    Dim
    for i in [10,20,40,80,160,320,640,900,1280,1700,2000,2560]:
        print i 
        DIMS.append(i)
        M = MATRICE[0:i,0:i]

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

    np.savetxt(home + workspace + folder + "__TimeReal"+str(9)+".tt", TimeReal)
    np.savetxt(home + workspace + folder + "__TimeComplex"+str(9)+".tt", TimeComplex)
    np.savetxt(home + workspace + folder + "__ErrReal"+str(9)+".tt", ErrReal)
    np.savetxt(home + workspace + folder + "__ErrComplex"+str(9)+".tt", ErrComplex)
    np.savetxt(home + workspace + folder + "__DIMS"+str(9)+".tt", DIMS)
