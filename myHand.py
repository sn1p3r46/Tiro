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

def findBlocks0(matrix):
    dim = matrix.shape[0]
    i=j=0; eb = sp.ndarray(dim,sp.integer); bb = sp.ndarray(dim,sp.integer)
    dims = sp.ndarray(dim,sp.int8); wb = sp.ndarray(dim)
    while (i<dim):
        if (i<dim-1 and matrix[i+1,i]!=0):
            bb[j] = i
            print str(bb[j])+"bb"
            eb[j] = i+1
            print str(eb[j])+"eb"
            dims[j] = 2; wb[i] = wb[i+1]=j
            i+=1
        else:
            bb[j] = eb[j] = i
            dims[j] = 1
            wb[i] = j
        i+=1
        j+=1
    #print bb; print"\n",print eb; print"\n",print wb; print"\n" ,print dims; print"\n"
    return bb,eb,dims,wb,i,j #i,wb dims maybe are useless


def getBlocks0(matrix):
    #add findBlocks0() return variables here!
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    #a = sp.ndarray(j)
    b = []
    for i in range(0,j):
        b.append(matrix[bb[i]:eb[i]+1,bb[i]:eb[i]+1]) #sp.append may be high cost operation
    return b

#much more simple! It Works!!! needs only: matrix, j(block number) and
#ordered dims of the blocks

def getBlocks1(matrix):
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    b = []
    k=0
    for i in range (0,j):
        b.append(matrix[k:k+dims[i],k:k+dims[i]])
        k += dims[i]
        print k
        print dims[i]; print '<--Dims'

    return b


#get eigenvalues perche' prendo quelli con parte immaginaria > 0 ?

def gev(matrix):
    x = (matrix[0,0]+matrix[1,1])/2
    y = sp.sqrt(-sp.square(matrix[0,0]-matrix[1,1])-4*matrix[1,0]*matrix[0,1])/2
    return complex(x,y)


# nroot = r^n(cos(ang*n) + i sin(ang*n))

def csr0(complex_n):
    s = sp.sqrt(sp.square(complex_n.real)+sp.square(complex_n.imag)) # np.absolute()
    angle = np.angle(complex_n)
    return sp.sqrt(s)*(complex((sp.cos(angle*0.5)), (sp.sin(angle*0.5))))


def csr1(complex_n):
    CMP0 = (sp.sqrt(sp.square(complex_n.real)+sp.square(complex_n.imag)))/2 #np.absolute()/2
    CMP1 = (complex_n.real/2)
    return complex(sp.sqrt(CMP0+CMP1),sp.sign(complex_n.imag)*(sp.sqrt(CMP0-CMP1)))


# From the book

def csr2(complex_n):
    t = sp.sqrt((sp.absolute(complex_n.real)+(sp.sqrt(sp.square(complex_n.real)+sp.square(complex_n.imag))))/2)
    if(complex_n.real >= 0):
        return complex(t,(complex_n.imag/(2*t)))
    else:
        return complex((complex_n.imag/(2*t)),t)

# De Moivre
    # -- square roots of a complex number are 2; n-roots of a square number are n
    # -- s^(1/n)* [cos((ang/2)+(2pi*k)/n) +i sin((ang/2)+(2pi*k)/n)] with k = [0...n-1]
    # http://www-thphys.physics.ox.ac.uk/people/FrancescoHautmann/Cp4/sl_clx_11_4_cls.pdf

def csr3(complex_n):
    ang = sp.angle(complex_n) # sp.arctan(a.imag/a.real) why it does not work?!?!
    r = sp.sqrt(sp.square(complex_n.real)+sp.square(complex_n.imag))
    if (sp.sin(ang/2)>=0): #sin>0
        return sp.sqrt(r)*(complex(sp.cos(ang/2),sp.sin(ang/2)))
    else:
        return sp.sqrt(r)*(complex(sp.cos((ang/2)+sp.pi),sp.sin((ang/2)+sp.pi)))

    #r1 = sp.sqrt(r)*(complex(sp.cos(ang/2),sp.sin(ang/2)))
    #r2 = sp.sqrt(r)*(complex(sp.cos((ang/2)+sp.pi),sp.sin((ang/2)+sp.pi)))
    #return r1,r2

def twobytworoot(matrix):

    if (matrix.shape[0]==2):
        a = csr3(gev(matrix)).real
        ris=sp.ndarray(shape=(2,2))
        ris[0,0] = a + (1/(4*a))*(matrix[0,0] - matrix[1,1])
        ris[1,1] = a - (1/(4*a))*(matrix[0,0] - matrix[1,1])
        ris[0,1] = (1/(2*a))*matrix[0,1]
        ris[1,0] = (1/(2*a))*matrix[1,0]
    else:
        print("sp.sqrt()\n")
        print(sp.sqrt(matrix))
        return sp.sqrt(matrix)
    return ris

def diagRoots(matrix):
    ris = sp.zeros(shape=matrix.shape,dtype=complex)
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    k=0
    for i in range (0,j):
        ris[k:k+dims[i],k:k+dims[i]] = twobytworoot(matrix[k:k+dims[i],k:k+dims[i]])
        print(matrix[k:k+dims[i],k:k+dims[i]])
        k += dims[i]
        print k
        print dims[i]; print '<--Dims'
    return ris

################################################################################
"""
def innerRoots(matrix):
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    for k in range():
"""        


################################################################################




if __name__ == "__main__":
    h = grism(10)
    """
    T, Z = applySchur(h)
    pmof(T,"BlockedMatrix.out")
    #applySchur(gmff("Source.out"))
    """

    print "\nImported\n"
    print "###########################"
