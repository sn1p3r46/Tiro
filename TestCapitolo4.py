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
    plt.plot(range(0,DeltaTime.shape[1]),(DeltaTime[5])-0.07,color="red",label="Matr dim: 320",lw=2)
    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[4],color="blue",label="Matr dim: 160",lw=2)
    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[3],color="black",label="Matr dim: 80",lw=2)
    plt.title("Block Dimensions Test on Small Matrices")
    plt.ylabel("Time")
    plt.xlabel("Blocks Dimension")
    plt.legend(loc='upper left')

    #plt.show()
    plt.show()
    """
    indici = DeltaTime.argmin(axis=1)
    X = range(0,DeltaTime.shape[0])

    plt.plot(X,indici,color="black",label="Dimensione dei blocchi ottimale rispetto alle matrici",lw=2)
    plt.ylabel("Time")
    plt.xlabel("Blocks Dimension")
    plt.show()

    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[7],color="black",label="1280",lw=2)
    plt.xlabel("Blocks Dimension")
    plt.ylabel("Time")
    plt.legend(loc='upper left')
    plt.show()
    """
    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[8]/270.0,color="red",label="Matr dim: 2560",lw=2)
    plt.xlabel("Blocks Dimension")
    plt.ylabel("Time")
    plt.legend(loc='upper left')
    #plt.show()

    plt.plot(range(0,DeltaTime.shape[1]),DeltaTime[6],color="black",label="Matr dim: 640",lw=2)
    plt.xlabel("Blocks Dimension")
    plt.ylabel("Time")
    plt.title("Block Dimensions Test on Big Matrices")
    plt.legend(loc='upper left')
    plt.show()

def plotfinal():
    ####################################################################################
    ##########################      CAPITOLO 4      ####################################
    ####################################################################################


    ##########################################
    ##############      TIME TEST
    ##########################################
    DeltatimeReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__TimeReal9.tt")
    DeltatimeComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__TimeComplex9.tt")
    MatrixDims = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__DIMS9.tt")
    X = range(0,DeltatimeReal.shape[0])
    #plt.xscale('log')
    plt.xlabel("Matrix Dimension")
    plt.ylabel("Time Sec.")
    plt.plot(MatrixDims,DeltatimeReal, color="red",label="Proposed Real Schur", lw=2)
    plt.plot(MatrixDims,DeltatimeComplex, color="blue",label="SciPy Complex Schur", lw=2)
    plt.title("Time Test Proposed vs SciPy")
    plt.legend(loc='upper left')
    plt.show()


    ##########################################
    ##############      ERROR TEST
    ##########################################

    plt.xlabel("Matrix Dimension")
    plt.ylabel("Relative Error")
    ErrReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__ErrReal9.tt")
    ErrComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__ErrComplex9.tt")
    plt.plot(MatrixDims,ErrReal, color="red",label="Proposed Real Schur", lw=2)
    plt.plot(MatrixDims,ErrComplex, color="blue",label="SciPy Complex Schur", lw=2)
    plt.title("Relative Error Test Proposed vs SciPy")
    plt.legend(loc='upper left')
    plt.show()

    plt.xlabel("Matrix Dimension")
    plt.ylabel("Time Sec.")
    DeltatimeReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/TestCapitolo2/Final9/__TimeReal9.tt")
    plt.plot(MatrixDims,DeltatimeReal, color="black",label="Real Schur", lw=2)
    DeltatimeComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__TimeComplex9.tt")
    plt.plot(MatrixDims,DeltatimeComplex, color="grey",label="SciPy Complex Schur", lw=2)
    plt.title("Relative Error Complex Schur (Chap2) vs SciPy")
    plt.title("Time Test")
    plt.legend(loc='upper left')
    plt.show()


    ####################################################################################
    ###########       CAPITOLO 2    ####################################################
    ####################################################################################


    ##########################################
    ##############      ERROR TEST
    ##########################################

    plt.xlabel("Matrix Dimension")
    plt.ylabel("Relative Error")
    ErrReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/TestCapitolo2/Final9/__ErrReal9.tt")
    ErrComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/TestCapitolo2/Final9/__ErrComplex9.tt")
    plt.plot(MatrixDims,ErrReal, color="red",label="Real Schur", lw=2)
    plt.plot(MatrixDims,ErrComplex, color="blue",label="Complex Schur", lw=2)
    plt.title("Relative Error Tests Real vs Complex Schur")
    plt.legend(loc='upper left')
    plt.show()

    print ErrComplex
    print ErrComplex

    ##########################################
    ##############      TIME TEST
    ##########################################

    plt.xlabel("Matrix Dimension")
    plt.ylabel("Time Sec.")

    DeltatimeReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/TestCapitolo2/Final9/__TimeReal9.tt")
    DeltatimeComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/TestCapitolo2/Final9/__TimeComplex9.tt")

    plt.plot(MatrixDims,DeltatimeReal, color="red",label="Real Schur", lw=2)
    plt.plot(MatrixDims,DeltatimeComplex, color="blue",label="Complex Schur", lw=2)
    plt.title("Time Tests Real vs Schur")

    print DeltatimeReal.shape
    print DeltatimeComplex.shape

    plt.legend(loc='upper left')
    plt.show()

    def matriceacaso(dimensione):
        return scipy.random.rand(dimensione,dimensione)
        T,Z = sp.linalg.schur(scipy.random.rand(dimensione,dimensione))

        aggiustala(T)
        return T


    ###########################################################################################################
    ####################################    TUTTI INSIEME  ####################################################
    ###########################################################################################################

    ##########################################
    ##############      TIME TEST
    ##########################################


    plt.xlabel("Matrix Dimension")
    plt.ylabel("Time Sec.")

    DeltatimeReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/TestCapitolo2/Final9/__TimeReal9.tt")
    DeltatimeComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/TestCapitolo2/Final9/__TimeComplex9.tt")
    plt.plot(MatrixDims,DeltatimeReal, color="red",label="Real Schur", lw=2)
    plt.plot(MatrixDims,DeltatimeComplex, color="blue",label="Complex Schur", lw=2)

    DeltatimeReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__TimeReal9.tt")
    DeltatimeComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__TimeComplex9.tt")
    plt.plot(MatrixDims,DeltatimeReal, color="black",label="Proposed Real Schur", lw=2)
    plt.plot(MatrixDims,DeltatimeComplex, color="green",label="SciPy Complex Schur", lw=2)
    plt.title("Time Test")
    plt.legend(loc='upper left')
    plt.show()


    ##########################################
    ##############      ERROR TEST
    ##########################################
    plt.xlabel("Matrix Dimension")
    plt.ylabel("Relative Error")
    ErrReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__ErrReal9.tt")
    ErrComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/ttest/Final9/__ErrComplex9.tt")

    plt.plot(MatrixDims,ErrReal, color="black",label="Proposed Real Schur", lw=2)
    plt.plot(MatrixDims,ErrComplex, color="green",label="SciPy Complex Schur", lw=2)

    ErrReal = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/TestCapitolo2/Final9/__ErrReal9.tt")
    ErrComplex = sp.loadtxt("/home/GFZNAS01/homes/galloni/Tiro/IOFiles/TestCapitolo2/Final9/__ErrComplex9.tt")
    plt.plot(MatrixDims,ErrReal, color="red",label="Real Schur", lw=2)
    plt.plot(MatrixDims,ErrComplex, color="blue",label="Complex Schur", lw=2)
    plt.title("Relative Error Test")
    plt.legend(loc='upper left')
    plt.show()







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
