import scipy as sp
import numpy as np
from handlerFile import *
from handlerBlocchi import *
from roots import *


def innerRoots(matrix):
    ris = diagRoots(matrix)
    global bb,eb,dims,wb,i,j 
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    for z in range(1,j):
        for t in range(0,j-z):
            rij = matrix[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1]
            tempsum = summ(ris,dims[t],dims[t+z],t,z)
            if (dims[t]==1 and dims[z+t]==1):
                tnoto = rij - tempsum
                ris[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1] = tnoto/(selectBlock(t,t,ris)+selectBlock(z+t,z+t,ris))

            elif (dims[t]==2 and dims[t+z]==1):
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
                tnoto = rij - tempsum
                ujj = selectBlock(z+t,z+t,ris)
                uii = selectBlock(t,t,ris)
                sysris = sp.ndarray((2,1))
                coeff = sp.zeros((2,2))
                coeff[0,0] = uii[0,0] + ujj[0,0]
                coeff[1,1] = uii[0,0] + ujj[1,1]
                coeff[0,1] = ujj[1,0]
                coeff[1,0] = ujj[0,1]
                sysris = np.linalg.solve(coeff,tnoto.transpose())
                print tempsum
                ris[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1] = sysris.transpose()

            else:
                tnoto = rij - tempsum
                ujj = selectBlock(z+t,z+t,ris)
                uii = selectBlock(t,t,ris)
                coeff = sp.ndarray((4,4))
                coeff[0,0] = uii[0,0]+ujj[0,0]; coeff[0,1] = uii[0,1]; coeff[0,2] = ujj[1,0]; coeff[0,3] = 0
                coeff[1,0] = uii[1,0]; coeff[1,1] = uii[1,1]+ujj[0,0]; coeff[1,2] = 0; coeff[1,3] = ujj[1,0]
                coeff[2,0] = ujj[0,1]; coeff[2,1] = 0; coeff[2,2] = uii[0,0]+ujj[1,1]; coeff[2,3] = uii[0,1]
                coeff[3,0] = 0; coeff[3,1] = ujj[0,1]; coeff[3,2]= uii[1,0]; coeff[3,3]=uii[1,1] + ujj[1,1]
                b = rij - tempsum; cc = sp.ndarray((4,1)); cc[0] = b[0,0]; cc[1] = b[1,0]; cc[2]=b[0,1]; cc[3]=b[1,1]

                temp = sp.zeros((2,2))
                FF = np.linalg.solve(coeff,cc)
                for k in range(0,2):
                    for h in range(0,2):
                        temp[h,k] = FF[2*k+h,0];
    #                      temp[0,0] = FF[0,0]; temp[0,1] = FF[2,0]; temp[1,0] = FF[1,0];temp[1,1] = FF[3,0]

                ris[bb[t]:eb[t]+1,bb[t+z]:eb[t+z]+1] = temp
                #ris[bb[t]:eb[t],bb[t+z]:eb[t+z]+1]
    return ris


def summ(matrix,d1,d2,t,z):
        if (d1==1 and d2==1):
            tempsum = 0
        elif(d1==2 and d2==1):
            tempsum = sp.zeros((2,1))
        elif(d1==1 and d2==2):
            tempsum = sp.zeros((1,2))
        else:
            tempsum = sp.zeros((2,2))
        for k in range(t+1,z+t):
            tempsum += matrix[bb[t]:eb[t]+1,bb[k]:eb[k]+1].dot(matrix[bb[k]:eb[k]+1,bb[z+t]:eb[z+t]+1])

        return tempsum
