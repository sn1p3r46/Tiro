import scipy as sp
import numpy as np
import time
from os.path import expanduser
import os

workspace = "/Tiro/IOFiles/"
home = expanduser("~") #home variable contains the user's home folder absolute path

def createWorkSpace():
    if( not os.path.isdir(home+workspace)):
        print "Path does not exists"
        os.makedirs(home+workspace)
        print 'Don\'t warry I\' ve just created it'
    else:
        print "File Folder is OK!"


#Get Matrix From File
def gmff(filename,path=None):
    if (path is None):
        return np.loadtxt(home+workspace+filename)
    elif(path==0):
        return np.loadtxt(filename)
        print filenamef
    else:
        if(path[len(path)-1]!="/"):
            path+="/"
        return np.loadtxt(path+filename)
    return 0


#Print Matrix on File
def pmof(Matrix,filename=None,path=None):
    if (filename is None):
        filename = time.strftime("%y.%m.%d-%H.%M.%S")+".out"
    if (path is None):
        path = home+workspace
    #np.savetxt(path+filename,Matrix,fmt='%8.3f') #%17.10f
    np.savetxt(path+filename,Matrix,fmt='%4.0f') #%17.10f
    return 0


#Generate Random Integer Square Matrix
def grism(Size):
    return np.random.randint(0, high=100, size=(Size,Size))


#Generate Random Real Square Matrix
def grrsm(Size):
    return np.random.uniform(-100000, high=100000, size=(Size,Size))


#############TO RUN ##############
createWorkSpace()
