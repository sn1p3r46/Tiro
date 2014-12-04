import numpy as np
import scipy as sp
from scipy import linalg
import os
from handlerFile import *
from os.path import expanduser


def main():
    print "\n###########################"


def applySchur(SourceMatrix):
    T, Z = sp.linalg.schur(SourceMatrix)
#   pmof(SourceMatrix,"Source.out")
    pmof(T, "Triangular.out")
    pmof(Z, "Unitary.out")
    return 0


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
    h=grrsm(7)
    #applySchur(gmff("Source.out"))

    print "\nImported\n"
