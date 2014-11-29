import numpy as np
import scipy as sp
from scipy import linalg
import os
from handlerFile import *
from os.path import expanduser

def main():
    print "\n###########################"

def applySchur(SourceMatrix):
    T,Z = sp.linalg.schur(SourceMatrix)
    #pmof(SourceMatrix,"Source.out")
    pmof(T,"Triangular.out")
    pmof(Z,"Unitary.out")
    return 0


if __name__ == "__main__":
    a=grrsm(7)
    applySchur(gmff("Source.out"))

    print "\nImported\n"
