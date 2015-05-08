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



def matriceacaso(dimensione):
    return scipy.random.rand(dimensione,dimensione)
    T,Z = sp.linalg.schur(scipy.random.rand(dimensione,dimensione))

    aggiustala(T)
    return T
