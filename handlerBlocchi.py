import scipy as sp

def findBlocks0(matrix):
    dim = matrix.shape[0]
    i=j=0; eb = sp.ndarray(dim,sp.integer); bb = sp.ndarray(dim,sp.integer)
    dims = sp.ndarray(dim,sp.int8); wb = sp.ndarray(dim)
    while (i<dim):
        if (i<dim-1 and matrix[i+1,i]!=0):
            bb[j] = i
            eb[j] = i+1
            dims[j] = 2; wb[i] = wb[i+1]=j
            i+=1

        else:
            bb[j] = eb[j] = i
            dims[j] = 1
            wb[i] = j

        i+=1
        j+=1

    return bb,eb,dims,wb,i,j #i,wb dims maybe are useless

def getBlocks0(matrix):
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    b = []
    for i in range(0,j):
        b.append(matrix[bb[i]:eb[i]+1,bb[i]:eb[i]+1])

    return b


def getBlocks1(matrix):
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    b = []
    k=0
    for i in range (0,j):
        b.append(matrix[k:k+dims[i],k:k+dims[i]])
        k += dims[i]

    return b


def selectBlock(m,n,matrix):
    bb,eb,dims,wb,i,j = findBlocks0(matrix)
    #print matrix[bb[m]:eb[m]+1,bb[n]:eb[n]+1]
    return matrix[bb[m]:bb[m]+dims[m],bb[n]:bb[n]+dims[n]]
