import numpy as np
import scipy as sp
from scipy import linalg
import os
from handlerFile import *
from os.path import expanduser



def applySchur(SourceMatrix):
    T, Z = sp.linalg.schur(SourceMatrix)
    #pmof(SourceMatrix,"Source.out")
    #pmof(T, "Triangular.out")
    #pmof(Z, "Unitary.out")
    return T,Z

def getBlocks0(matrix):
    dim = matrix.shape[0]
    i=j=0; bb = eb = sp.ndarray(dim,sp.integer)
    dims = sp.ndarray(dim,sp.int8); wb = sp.ndarray(dim)
    while (i<dim):
        if (i<dim-1 and matrix[i+1,i]!=0):
            bb[j] = i
            print bb[j]
            eb[j] = i+1
            print eb[j]
            dims[j] = 2; wb[i] = wb[i+1]=j
            i+=1
        else:
            bb[j] = eb[j] = i
            dims[i] = 1
            wb[i] = j
        i+=1
        j+=1
    return bb,eb,dims,wb

def csr0(complex_n):
    s = sp.sqrt(sp.square(complex_n.real)+sp.square(complex_n.imag)) # np.absolute()
    angle = np.angle(complex_n)
    return sp.sqrt(s)*(complex((sp.cos(angle*0.5)), (sp.sin(angle*0.5))))


def csr1(complex_n):
    CMP0 = (sp.sqrt(sp.square(complex_n.real)+sp.square(complex_n.imag)))/2 #np.absolute()/2
    CMP1 = (complex_n.real/2)
    return complex(sp.sqrt(CMP0+CMP1),sp.sign(complex_n.imag)*(sp.sqrt(CMP0-CMP1)))


def csr2(complex_n):
    t = sp.sqrt((sp.absolute(complex_n.real)+(sp.sqrt(sp.square(complex_n.real)+sp.square(complex_n.imag))))/2)
    if(complex_n.real >= 0):
        return complex(t,(complex_n.imag/(2*t)))
    else:
        return complex((complex_n.imag/(2*t)),t)

if __name__ == "__main__":
    h = grism(10)
    """
    T, Z = applySchur(h)
    pmof(T,"BlockedMatrix.out")
    #applySchur(gmff("Source.out"))
    """

    print "\nImported\n"
    print "###########################"
