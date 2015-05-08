import numpy as np
import scipy as sp
import matplotlib
from matplotlib import pyplot as plt
import scipy.linalg as spl

TimeStart = np.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttestBlocks/TimeStart9.tt")
TimeEnd = np.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttestBlocks/TimeEnd9.tt")
DeltaTime = TimeEnd - TimeStart




########################---- TEST GRANDEZZA BLOCCHI ----##############################

def plotBlocchi():
    TimeStart = np.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttestBlocks/TimeStart9.tt")
    TimeEnd = np.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttestBlocks/TimeEnd9.tt")
    DeltaTime = TimeEnd - TimeStart

    X = DeltaTime.shape[1]
    #plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[8],color="green",label="8")
    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[4],color="brown",label="4")
    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[2],color="orange",label="6")
    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[3],color="black",label="5")

    plt.legend(loc='upper left')

    #plt.show()
    plt.show()

    indici = DeltaTime.argmin(axis=1)
    X = range(0,DeltaTime.shape[0])

    plt.plot(X,indici,color="black",label="Dimensione dei blocchi ottimale rispetto alle matrici")
    plt.ylabel("Time")
    plt.xlabel("Blocks Dimension")
    plt.show()

    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[7],color="black",label="1280")
    plt.xlabel("Blocks Dimension")
    plt.ylabel("Time")
    plt.legend(loc='upper left')
    plt.show()

    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[8],color="black",label="2560")
    plt.xlabel("Blocks Dimension")
    plt.ylabel("Time")
    plt.legend(loc='upper left')
    plt.show()

    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[6],color="black",label="640")
    plt.xlabel("Blocks Dimension")
    plt.ylabel("Time")
    plt.legend(loc='upper left')
    plt.show()

####################################################################################
####################################################################################

DeltatimeReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/9/DeltaTime_9.tt")
DeltatimeComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/9/ComplexDeltaTime1_9.tt")
MatrixDims = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/9/dimensionematricigenerate_9.tt")
X = range(0,DeltatimeReal.shape[0])

DeltatimeReal =  np.insert(DeltatimeReal,8,39.811653)
DeltatimeComplex =  np.insert(DeltatimeComplex,7,43.93453599999998)
MatrixDims =  np.insert(MatrixDims,7,900)

#plt.xscale('log')
plt.xlabel("Matrix Dimension")
plt.ylabel("Time")
plt.legend(loc='upper left')

plt.plot(MatrixDims,DeltatimeReal, color="black",label="Time Real")
plt.plot(MatrixDims,DeltatimeComplex, color="green",label="Time Complex")

plt.legend(loc='upper left')
plt.show()





def matriceacaso(dimensione):
    return scipy.random.rand(dimensione,dimensione)
    T,Z = sp.linalg.schur(scipy.random.rand(dimensione,dimensione))

    aggiustala(T)
    return T



"""

In [94]: FFF = MMM[0:640,0:640]

In [95]: timeTestRealSchur(FFF)

Out[95]: (14.858203, 1.8190959636500577e-16)

In [96]: timeTestComplexSchur(FFF)
Out[96]: (15.523201999999998, 2.1120068302034346e-16)

In [97]: FFF = MMM[0:900,0:900]

In [98]: timeTestRealSchur(FFF)
Out[98]: (39.811653, 1.8280833493050758e-16)

In [99]: timeTestComplexSchur(FFF)
Out[99]: (43.93453599999998, 2.2094605232873023e-16)

In [100]: FFF = MMM[0:1200,0:1200]

In [101]: timeTestRealSchur(FFF)
Out[101]: (95.017705, 1.8507274218634707e-16)

In [102]: timeTestComplexSchur(FFF)
Out[102]: (110.53159400000001, 2.3165688557369921e-16)

MMM = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/BigBlockMatrix5000.out"

"""
