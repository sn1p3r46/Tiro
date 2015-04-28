import numpy as np
import scipy as sp
from scipy import linalg
import os
from handlerFile import *
from os.path import expanduser
import pdb

def applySchur(SourceMatrix):
    T, Z = sp.linalg.schur(SourceMatrix)
    return T,Z

def findBlocks0(matrix):
    dim = matrix.shape[0]
    i=j=0; eb = sp.ndarray(dim,sp.integer); bb = sp.ndarray(dim,sp.integer)
    dims = sp.ndarray(dim,sp.int8); wb = sp.ndarray(dim)
    while (i<dim):
        if (i<dim-1 and matrix[i+1,i]!=0):
            bb[j] = i
            #print str(bb[j])+"bb"
            eb[j] = i+1
            #print str(eb[j])+"eb"
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
        #print(sp.sqrt(matrix))
        return sp.sqrt(matrix)
    return ris

def diagRoots(matrix):
    ris = sp.zeros(shape=matrix.shape)
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    k=0
    for i in range (0,j):
        print matrix[k:k+dims[i],k:k+dims[i]]
        ris[k:k+dims[i],k:k+dims[i]] = twobytworoot(matrix[k:k+dims[i],k:k+dims[i]])
        #print(matrix[k:k+dims[i],k:k+dims[i]])
        k += dims[i]
        #print k
        #print dims[i]; print '<--Dims'
    return ris

################################################################################

def innerRoots(matrix):
    ris = diagRoots(matrix)
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    for z in range(1,j):
        for t in range(0,j-z):
            rij = matrix[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1]
            #print rij
            #print (rij)
            if (dims[t]==1 and dims[z+t]==1):
                tempsum = 0
                tempsum0 = 0
                for k in range(t+1,z+t):
                    tempsum += ris[bb[t]:eb[t]+1,bb[k]:eb[k]+1].dot(ris[bb[k]:eb[k]+1,bb[z+t]:eb[z+t]+1])

                print("Tempsum(1x1):")
                print tempsum
                tnoto = rij - tempsum
                ris[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1] = tnoto/(selectBlock(t,t,ris)+selectBlock(z+t,z+t,ris))

            elif (dims[t]==2 and dims[t+z]==1):
                tempsum = sp.zeros((2,1))
                for k in range(t+1,t+z):
                    tempsum += ris[bb[t]:eb[t]+1,bb[k]:eb[k]+1].dot(ris[bb[k]:eb[k]+1,bb[z+t]:eb[z+t]+1])
                print("Tempsum(2x1):")
                print tempsum
                    #tempsum0 += selectBlock(t,k,ris).dot(selectBlock(k,z,ris))
                tnoto = rij - tempsum
                ujj = selectBlock(z+t,z+t,ris)
                uii = selectBlock(t,t,ris)
                sysris = sp.ndarray((2,1))
                coeff = sp.array(uii)
                coeff[0,0] += ujj[0,0]
                coeff[1,1] += ujj[0,0]
                sysris = np.linalg.solve(coeff,tnoto)
                ris[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1] = sysris

            elif (dims[t]==1 and dims[t+z]==2):
                tempsum = sp.zeros((1,2))
                for k in range(t+1,t+z):
                    tempsum += ris[bb[t]:eb[t]+1,bb[k]:eb[k]+1].dot(ris[bb[k]:eb[k]+1,bb[z+t]:eb[z+t]+1])
                #pdb.set_trace()
                print("Tempsum(1x2):")
                print tempsum
                tnoto = rij - tempsum
                ujj = selectBlock(z+t,z+t,ris)
                uii = selectBlock(t,t,ris)
                sysris = sp.ndarray((2,1))
                coeff = sp.zeros((2,2))
                print uii[0,0]
                coeff[0,0] = uii[0,0] + ujj[0,0]
                coeff[1,1] = uii[0,0] + ujj[1,1]
                coeff[0,1] = ujj[1,0]
                coeff[1,0] = ujj[0,1]
                sysris = np.linalg.solve(coeff,tnoto.transpose())
                print "tnoto.transpose()"
                print tnoto.transpose()
                print z
                print t
                print coeff
                print("SysRis")
                print sysris.shape
                print matrix[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1]
                print("SysRis")
                ris[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1] = sysris.transpose()
            else:
                print("sono entrato nell'else")
                tempsum = sp.zeros((2,2))
                for k in range(t+1,t+z):
                    tempsum += ris[bb[t]:eb[t]+1,bb[k]:eb[k]+1].dot(ris[bb[k]:eb[k]+1,bb[z+t]:eb[z+t]+1])
                tnoto = rij - tempsum
                print tnoto
                print("Tempsum(2x2):")
                print tempsum
                ujj = selectBlock(z+t,z+t,ris)
                uii = selectBlock(t,t,ris)
                coeff = sp.ndarray((4,4))
                print uii
                print ujj
                print coeff
                coeff[0,0] = uii[0,0]+ujj[0,0]; coeff[0,1] = uii[0,1]; coeff[0,2] = ujj[1,0]; coeff[0,3] = 0
                coeff[1,0] = uii[1,0]; coeff[1,1] = uii[1,1]+ujj[0,0]; coeff[1,2] = 0; coeff[1,3] = ujj[1,0]
                coeff[2,0] = ujj[0,1]; coeff[2,1] = 0; coeff[2,2] = uii[0,0]+ujj[1,1]; coeff[2,3] = uii[0,1]
                coeff[3,0] = 0; coeff[3,1] = ujj[0,1]; coeff[3,2]= uii[1,0]; coeff[3,3]=uii[1,1] + ujj[1,1]
                b = rij - tempsum; cc = sp.ndarray((4,1)); cc[0] = b[0,0]; cc[1] = b[1,0]; cc[2]=b[0,1]; cc[3]=b[1,1]
                print cc
                print coeff
                temp = sp.zeros((2,2))
                FF = sp.linalg.solve(coeff,cc)
                for k in range(0,2):
                    for h in range(0,2):
                        temp[h,k] = FF[2*k+h,0];
    #                          temp[0,0] = FF[1,0]; temp[0,1] = FF[2,0]; temp[1,0] = FF[1,0];temp[1,1] = FF[3,0]
                print FF
                print temp
                ris[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1] = temp
                #ris[bb[t]:eb[t],bb[t+z]:eb[t+z]+1]
    return ris


def innerRoots1(matrix):
    bb,eb,dims,wb,v,nBlocks = findBlocks0(matrix)
    ris = diagRoots(matrix)
    for j in range(1, nBlocks):
        for i in range (0, nBlocks - j):
            if(dims[i]==1 and dims[i+j]==1):
                tnoto = tnotocaso1(i,j,ris)
                rij = selectBlock(i,i+j,matrix)
                uii = selectBlock(i,i,ris)
                ujj = selectBlock(i+j,i+j,ris)
                ris[bb[i]:eb[i],bb[i+j]:eb[i+j]+1] = (rij - tnoto)/(uii*ujj)

            elif(dims[i]==2 and dims[i+j]==1):
                tnoto = tnotocaso2(i,j,ris,dims[i],dims[i+j])
                rij = selectBlock(i,i+j,matrix)
                uii = selectBlock(i,i,ris)
                ujj = selectBlock(i+j,i+j,ris)
                coeff = sp.ndarray((2,2),dtype=complex)
                coeff[0,0] = uii[0,0] + ujj[0,0]
                coeff[1,1] = uii[1,1] + ujj[0,0]
                coeff[1,0] = uii[1,0]
                coeff[0,1] = uii[0,1]
                print i,j, ris[bb[i]:eb[i+j],bb[i+j]:eb[i+j]+1]
                ris[bb[i]:eb[i],bb[i+j]:eb[i+j]+1] = sp.linalg.solve(coeff,rij-tnoto)

            elif(dims[i]==1 and dims[i+j]==2):
                tnoto = tnotocaso2(i,j,ris,dims[i],dims[i+j])
                rij = selectBlock(i,i+j,matrix)
                uii = selectBlock(i,i,ris)
                ujj = selectBlock(i+j,i+j,ris)
                coeff = sp.ndarray((2,2), dtype=complex)
                coeff[0,0] = uii[0,0] + ujj[0,0]
                coeff[1,1] = uii[0,0] + ujj[1,1]
                coeff[1,0] = ujj[0,1]
                coeff[0,1] = ujj[1,0]
                print "coeff"
                print coeff
                print (rij-tnoto).transpose()
                ris[bb[i]:eb[i],bb[i+j]:eb[i+j]+1] = sp.linalg.solve(coeff.transpose(),(rij-tnoto).transpose())

            elif(dims[i]==2 and dims[i+j]==2):
                tnoto = tnotocaso2(i,j,ris,dims[i],dims[i+j])
                rij = selectBlock(i,i+j,matrix)
                uii = selectBlock(i,i,ris)
                ujj = selectBlock(i+j,i+j,ris)
                coeff = sp.ndarray((4,4), dtype=complex)
                coeff[0,0] = uii[0,0]+ujj[0,0]; coeff[0,1] = uii[1,1]; coeff[0,2] = ujj[1,0]; coeff[0,3] = 0
                coeff[1,0] = uii[1,0]; coeff[1,1] = uii[1,1]+ujj[0,0]; coeff[1,2] = 0; coeff[1,3] = ujj[1,0]
                coeff[2,0] = ujj[0,1]; coeff[2,1] = 0; coeff[2,2] = uii[0,0]+ujj[1,1]; coeff[2,3] = uii[0,1]
                coeff[3,0] = 0; coeff[3,1] = ujj[0,1]; coeff[3,2]= uii[1,0]; coeff[3,3]=uii[1,1] + ujj[1,1]
                ###############################
                b = rij-tnoto; bb = sp.ndarray(4); bb[0] = b[0,0]; bb[1] = b[0,1]; bb[2]=b[1,0]; bb[3]=b[1,1]
                ris[bb[i]:eb[i],bb[i+j]:eb[i+j]+1] = sp.linalg.solve(coeff,bb)
            else:
                break

    return ris


def tnotocaso2(i,j,matrix,x,y):
    tnoto = sp.zeros((x,y), dtype=complex)
    for k in range(i+1,i+j):
        tnoto += selectBlock(i,k,matrix).dot(selectBlock(k,i+j,matrix))
    print tnoto
    return tnoto


def tnotocaso1(i,j,matrix):
    tnoto = 0
    for k in range(i+1,i+j):
        tnoto += selectBlock(i,k,matrix).dot(selectBlock(k,i+j,matrix))
    return tnoto

        #print "<------------------------/n"


def selectBlock(m,n,matrix):
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    #print matrix[bb[m]:eb[m]+1,bb[n]:eb[n]+1]
    return matrix[bb[m]:bb[m]+dims[m],bb[n]:bb[n]+dims[n]]


################################################################################


if __name__ == "__main__":
    # h = grism(10)
#    T, Z = applySchur(h)
#    pmof(T,"BlockedMatrix.out")
    #applySchur(gmff("Source.out"))
    print "\nImported\n"
    print "###########################"
