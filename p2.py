import myHand

def gcm():
    while(True):
        a = grism(7)
        A,B = applySchur(a)
        if (getblocks1(A)!=False):
            return A

def getblocks1(matrix):
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
              if(matrix[i,i]<0):
                  return False
          i+=1
          j+=1
      return bb,eb,dims,wb,i,j
    #i wb dims maybe are useless
