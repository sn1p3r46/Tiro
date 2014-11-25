import scipy as sp
import numpy as np
from os.path import expanduser
import os

workspace = "/Tiro/IOFiles/"
home = expanduser("~") #home variable contains the home folder's absolute path

def createWorkSpace():
    if( not os.path.isdir(home+workspace)):
        print "Path does not exists"
        os.makedirs(home+workspace)
        print 'Don\'t warry I\' ve just created it'
    else:
        print "File Folder is OK!"


#Give Me Matrix From File
def gmmff():
    return 0


#Print Matrix on File
def pmof(Matrix,filename=None):
    if (filename is None):
        #s
        filename = time.strftime("%y.%m.%d-%H.%M.%S")+".out"
    np.savetxt(home+workspace+filename,Matrix,fmt='%17.10f')
    return 0


#Generate Random Integer Square Matrix
def grism(Size):
    return np.random.randint(-1000, high=1000, size=(Size,Size))


#Generate Random Real Square Matrix
def grrsm():
    return 0


#############TO RUN ##############
createWorkSpace()
