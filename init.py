
from os.path import expanduser
import os

import scipy as sp
import scipy.linalg

import time

workspace = "/Tiro/IOFiles/"
home = expanduser("~")

ttest = "/ttest/"
ttestBlocks = "/ttestBlocks/"
ttest2 = "/ttest2/"
execfile("myHand.py")
execfile("sqrtm5.py")
execfile("sqrtm4.py")
execfile("aggiustala.py")
A = gmff("BlockedMatrix.out")

"""
def block_structure(T):
    n = len(T)
    tol = 1e-15
    i,j = 0,0
    bb = sp.empty(n,dtype="int")
    eb = sp.empty(n,dtype="int")
    s  = sp.empty(n,dtype="int")
    while i < n-1:
        bb[j] = i
        if abs(T[i+1,i])<tol:
            i +=1
            s[j] = 1
            eb[j] = i
        else:
            i +=2
            s[j] = 2
            eb[j] = i
        j += 1
    if i == n-1:
        bb[j],eb[j] = i,i+1
        j+=1
    return j, bb, eb, s
"""

def timeTest(ii):
    I = np.zeros(ii)
    F = np.zeros(ii)
    dims = []
    Err = []
    I1 = np.zeros(ii)
    F1 = np.zeros(ii)
    dims1 = []
    Err1 = []
    #T,Z = sp.linalg.schur(B)
    #aggiustala(T)
    Init_Dim = 10
    for i in range(0,ii):
        B = scipy.random.randint(-100,100,(Init_Dim,Init_Dim))
        print "\nCreo matrice B di dimensione: " + str(Init_Dim)+"\n"
        print "\nOra calcolo schur di B!"
        T,Z = sp.linalg.schur(B)
        print "\nDecomposizione di schur terminata con sucesso !!!\n "
        aggiustala(T)
        print "\nAggiustata ! ! ! Ora calcolo sqrtm5(T)\n"
        I[i]=time.clock()
        MM = sqrtm5(T)
        F[i]=time.clock()
        print "\nCalcolata! Ora calcolo sp.linalg.sqrtm(T)\n"
        Err.append(sp.linalg.norm(MM.dot(MM)-T)/sp.linalg.norm(T))
        dims.append(Init_Dim)

        I1[i]=time.clock()
        MM = sp.linalg.sqrtm(T)
        F1[i]=time.clock()
        Err1.append(sp.linalg.norm(MM.dot(MM)-T)/sp.linalg.norm(T))
        dims1.append(Init_Dim)
        print "\nFinita dimensione--------------------" + str(Init_Dim) +"\n"
        Init_Dim=Init_Dim*2



    np.savetxt(home + workspace + ttest + "TimeStart_"+str(ii)+".tt", I, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' inizio del calcolo per sqrtm5()')
    np.savetxt(home + workspace + ttest + "TimeEnd_"+str(ii)+".tt", F, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' fine del calcolo per sqrtm5()')
    np.savetxt(home + workspace + ttest + "DeltaTime_"+str(ii)+".tt", (F-I), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Tempo realmente impiegato da sqrtm5,\n il numero ii indica la dim della matrice come 10*2^ii con ii che va da 0 a ii')
    np.savetxt(home + workspace + ttest + "dimensionematricigenerate_"+str(ii)+".tt", dims, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' dim delle matrici su cui sono effettuati i calcoli')
    np.savetxt(home + workspace + ttest + "ErrRelativo_"+str(ii)+".tt", Err, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Errore Relativo Per ogni matrice')


    np.savetxt(home + workspace + ttest + "TimeStart1_"+str(ii)+".tt", I1, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' inizio del calcolo per sp.linalg.sqrtm()')
    np.savetxt(home + workspace + ttest + "TimeEnd1_"+str(ii)+".tt", F1, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' fine del calcolo per sp.linalg.sqrtm()')
    np.savetxt(home + workspace + ttest + "DeltaTime1_"+str(ii)+".tt", (F1-I1), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Tempo realmente impiegato da sp.linalg.sqrtm(),\n il numero ii indica la dim della matrice come 10*2^ii con ii che va da 0 a ii')
    np.savetxt(home + workspace + ttest + "dimensionematricigenerate1_"+str(ii)+".tt", dims1, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' dim delle matrici su cui sono effettuati i calcoli')
    np.savetxt(home + workspace + ttest + "ErrRelativo1_"+str(ii)+".tt", Err1, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Errore Relativo Per ogni matrice')


    print "I:"
    print I
    print "F:"
    print F
    print "F - I:"
    print F - I
    print "Dims:"
    print dims
    print "I1:"
    print I1
    print "F1:"
    print F1
    print "F1 - I1:"
    print F1 - I1
    print "Dims:"
    print dims1
    print "ASSOLUTO"
    print Err
    print Err1

def timeTest1(ii,jj):
    I = np.zeros((ii,jj-2))
    F = np.zeros((ii,jj-2))
    print I.shape
    dims = []
    #T,Z = sp.linalg.schur(B)
    #aggiustala(T)
    Init_Dim = 10
    for i in range(0,ii):
        B = scipy.random.randint(-100,100,(Init_Dim,Init_Dim))
        T,Z = sp.linalg.schur(B)
        aggiustala(T)
        for j in range(2,jj):
            I[i,j-2] = time.clock()
            MM   = sqrtm5(T,j)
            F[i,j-2] = time.clock()
	dims.append(Init_Dim)
	Init_Dim=Init_Dim*2

    FI = F - I

# il seguente output e' incentrato a verificare se ci sono cambiamenti a livello di tempo al variare delle dimensioni dei blocchi grandi insiame alle matrici
# ogni colonna contiene il tempo al variare dei blocchi sulla matrice,
# mentre ogni riga aumenta la dimensione della matrice in ordine quadratico;
# min e max time  invece misurano il massimo cambiamento di tempo per ogni riga, per vedere con deltatime poi se effettivamente la dimensione dei blocchi conta ai fini della velocita.

    np.savetxt(home + workspace + ttestBlocks + "TimeStart"+str(ii)+".tt", I, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' inizio del calcolo per sqrtm5()')
    np.savetxt(home + workspace + ttestBlocks + "TimeEnd"+str(ii)+".tt", F, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' fine del calcolo per sqrtm5()')
    np.savetxt(home + workspace + ttestBlocks + "DeltaTime"+str(ii)+".tt", (F-I), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Tempo realmente impiegato da sqrtm5,\n il numero ii indica la dim della matrice come 10*2^ii con ii che va da 0 a ii')
    np.savetxt(home + workspace + ttestBlocks + "dimensionematricigenerate"+str(ii)+".tt", dims, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' dim delle matrici su cui sono effettuati i calcoli')
#   np.savetxt(home + workspace + ttestBlocks + "ErrRelativo"+str(ii)+".tt", Err, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Errore Relativo Per ogni matrice')
    np.savetxt(home + workspace + ttestBlocks + "MaxTime" +str(ii)+".tt", np.amax(FI, axis=1), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' tempo massimo al variare delle dims dei blocchi per ogni matrice')
    np.savetxt(home + workspace + ttestBlocks + "MinTime" +str(ii)+".tt", np.amin(FI, axis=1), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' tempo minimo al variare delle dimens dei blocchi per ogni matrice')
    np.savetxt(home + workspace + ttestBlocks + "DeltaTime"+ str(ii)+".tt", (np.amax(FI, axis=1) - np.amin(FI, axis=1)) , fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' differenza tra il tempo massimo ed il tempo minimo per ogni matrice')

    print "I:"
    print I
    print "F:"
    print F
    print "F-I:"
    print F - I
    print "dims:"
    print dims
    print "MAX:"
    print np.amax(FI, axis=1)
    print "MIN:"
    print np.amin(FI, axis=1)
    print "DIFF:"
    print np.amax(FI, axis=1) - np.amin(FI, axis=1)


def timeTest2(ii=1,dimensione=10240):
    I = np.zeros(ii)
    F = np.zeros(ii)
    dims = []
    Err = []
    I1 = np.zeros(ii)
    F1 = np.zeros(ii)
    dims1 = []
    Err1 = []
    #T,Z = sp.linalg.schur(B)
    #aggiustala(T)

    for i in range(ii-1,ii):
        B = scipy.random.rand(dimensione,dimensione)
        print "\nCreo matrice B di dimensione: " + str(dimensione)+"\n"
	print "\nOra calcolo schur di B!\n"
        T,Z = sp.linalg.schur(B)
	del B
	print "\nDecomposizione di schur terminata con sucesso !!! AGGIUSTO\n "
        aggiustala(T)
        print "\nAggiustata ! ! ! Ora calcolo sqrtm5(T)\n"
        I[i]=time.clock()
        MM = sqrtm5(T)
        F[i]=time.clock()
        print "\nCalcolata! Ora calcolo sp.linalg.sqrtm(T)\n"
        Err.append(sp.linalg.norm(MM.dot(MM)-T)/sp.linalg.norm(T))
        dims = [dimensione]
        I1[i]=time.clock()
        MM = sp.linalg.sqrtm(T)
        F1[i]=time.clock()
        Err1.append(sp.linalg.norm(MM.dot(MM)-T)/sp.linalg.norm(T))
        dims1 = [dimensione]
	print "\nFinita dimensione--------------------" + str(dimensione) +"\n"



        np.savetxt(home + workspace + ttest2 + "00TimeStart_"+str(ii)+".tt", I, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' inizio del calcolo per sqrtm5()')
        np.savetxt(home + workspace + ttest2 + "00TimeEnd_"+str(ii)+".tt", F, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' fine del calcolo per sqrtm5()')
        np.savetxt(home + workspace + ttest2 + "00DeltaTime_"+str(ii)+".tt", (F-I), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Tempo realmente impiegato da sqrtm5,\n il numero ii indica la dim della matrice come 10*2^ii con ii che va da 0 a ii')
        np.savetxt(home + workspace + ttest2 + "00dimensionematricigenerate_"+str(ii)+".tt", dims, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' dim delle matrici su cui sono effettuati i calcoli')
        np.savetxt(home + workspace + ttest2 + "00ErrASSOLUTO"+str(ii)+".tt", Err, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Errore Relativo Per ogni matrice')


        np.savetxt(home + workspace + ttest2 + "00TimeStart1_"+str(ii)+".tt", I1, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' inizio del calcolo per sp.linalg.sqrtm()')
        np.savetxt(home + workspace + ttest2 + "00TimeEnd1_"+str(ii)+".tt", F1, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' fine del calcolo per sp.linalg.sqrtm()')
        np.savetxt(home + workspace + ttest2 + "00DeltaTime1_"+str(ii)+".tt", (F1-I1), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Tempo realmente impiegato da sp.linalg.sqrtm(),\n il numero ii indica la dim della matrice come 10*2^ii con ii che va da 0 a ii')
        np.savetxt(home + workspace + ttest2 + "00dimensionematricigenerate1_"+str(ii)+".tt", dims1, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' dim delle matrici su cui sono effettuati i calcoli')
        np.savetxt(home + workspace + ttest2 + "00ErrASSOLUTO_"+str(ii)+".tt", Err1, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments=' Errore Relativo Per ogni matrice')


        print "I:"
        print I
        print "F:"
        print F
        print "F - I:"
        print F - I
        print "Dims:"
        print dims
        print "I1:"
        print I1
        print "F1:"
        print F1
        print "F1 - I1:"
        print F1 - I1
        print "Dims:"
        print dims1
        print "ASSOLUTO"
        print Err
        print Err1

def matriceacaso(dimensione):
    return scipy.random.rand(dimensione,dimensione)
    T,Z = sp.linalg.schur(scipy.random.rand(dimensione,dimensione))
    
    aggiustala(T)
    return T
#print min(sp.linalg.eigvals(T))
