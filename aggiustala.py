def aggiustala(M):
    m,bb,eb,s = block_structure(M)
    for i in range(0,m):
        X = M[bb[i]:eb[i],bb[i]:eb[i]]
        if(X.shape[0]==1):
            M[bb[i]:eb[i],bb[i]:eb[i]] = abs(M[bb[i]:eb[i],bb[i]:eb[i]])
        else:
            X[0,0] = abs(X[0,0])
            X[1,1] = abs(X[0,0])
            X[0,1] = -abs(X[0,1])
            X[1,0] = abs(X[0,1])
