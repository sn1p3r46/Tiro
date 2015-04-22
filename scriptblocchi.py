

import scipy as sp

import scipy.linalg

from scipy.linalg.lapack import ztrsyl, dtrsyl

import scipy.linalg as spl



D = sp.random.randint(-100,100,(2,8000))


def blocks1(M):
    for i in range(0,2000,+4):
        Uii = M[:,i:i+2]
        Ujj = M[:,i+2:i+4]

        Tnoto = sp.random.randint(-100,100,(2,2))

        #shapeUii = Uii.shape[0]
        #shapeUjj = Ujj.shape[0]

        x, scale, info = dtrsyl(Uii, Ujj, Tnoto)
        ris = x*scale

def blocks2(M):
    for i in range(0,2000,+4):
        Uii = M[:,i:i+2]
        Ujj = M[:,i+2:i+4]
        Tnoto = sp.random.randint(-100,100,(2,2))

        shapeUii = Uii.shape[0]
        shapeUjj = Ujj.shape[0]

        Tnoto = Tnoto.reshape((shapeUii*shapeUjj),1,order="F")
        ris =   sp.eye(shapeUii)
        ris2=sp.eye(shapeUjj)
        #        sp.linalg.solve(sp.kron(Ujj,Uii) +
        #        sp.kron(Ujj,Uii),Tnoto).reshape(shapeUii,shapeUjj,order="F")
        #        sp.linalg.solve(sp.kron(sp.eye(shapeUjj),Uii) +
        #        sp.kron(Ujj.T,sp.eye(shapeUii)),
        # Tnoto).reshape(shapeUii,shapeUjj,order="F")




def INsc1(A):
    X = np.copy(A)
    E = 0.5*(np.eye(A.shape[0])-A)
    dA = spl.det(A)**(0.5)
    for i in range (1,5):
        print i
        detX =  spl.det(X)
        uk = np.abs(detX/dA)**(-1.0/i)
        Etk = (E + (0.5*X))/uk - (0.5)*uk*X
        X = uk*X + Etk
        print X
        E = (-0.5)*(Etk.dot(spl.inv(X))).dot(Etk)
    return X
