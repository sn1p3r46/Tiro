
import scipy as sp
import scipy.linalg

execfile("myHand.py")
execfile("sqrtm5.py")
execfile("sqrtm4.py")
execfile("aggiustala.py")
A = gmff("BlockedMatrix.out")

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

B = scipy.random.randint(-100,100,(1200,1200))

T,Z = sp.linalg.schur(B)

aggiustala(T)
print min(sp.linalg.eigvals(T))
