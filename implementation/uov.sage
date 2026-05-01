def Upper(M):
    '''
    returns upper triangular matrix Upper(M), defined as in MAYO specification
    '''
    n = M.nrows()
    FF = M.base_ring()
    upper = matrix(FF, n)
    for i in range(n):
        upper[i, i] = M[i, i]
    for i in range(n):
        for j in range(i+1, n):
            upper[i, j] = M[i, j] + M[j, i]    
    return upper

def UpperRandom(FF, n): 
    '''
    returns random upper triangular matrix M
    '''
    M = matrix(FF, n)
    for i in range(n):
        for j in range(i, n):
            M[i, j] = FF.random_element()  
    return M

def KeyGen(v, o, m, q): 
    '''
    returns key pair of UOV
    this code is for research purposes only and it is not secure for real-world cryptographic applications
    '''
    n = o + v    
    FF = GF(q)

    I1 = identity_matrix(FF, v)
    I2 = identity_matrix(FF, o)
    O = matrix(FF, o, v)
    B = random_matrix(FF, v, o)
    S = block_matrix([[I1, -B], [O, I2]]) # compact form of S, Oil subspace = span [B I2]

    P = []
    for i in range(m):
        
        if FF.characteristic() != 2 : # generate key pair in a way similar of QR-UOV (no quotient ring structure here)
            Fi = matrix(FF, n)
            while Fi.determinant() == 0: 
                F0 = matrix(FF, o)             
                F1 = random_matrix(FF, v)
                F1 = F1 + F1.transpose()
                F2 = random_matrix(FF, v, o)
                Fi = block_matrix([[F1, F2], [F2.transpose(), F0]])
            Pi = S.transpose()*Fi*S
            P.append(Pi)
            
        else: # generate key pair of UOV/MAYO
            Pi = matrix(FF, n)
            while Pi.determinant() == 0: 
                P0 = matrix(FF, o, v)
                P1 = UpperRandom(FF, v)
                P2 = random_matrix(FF, v, o)
                P3 = Upper(-B.transpose()*P1*B + B.transpose()*P2)
                Pi = block_matrix([[P1, P2], [P0, P3]])
            P.append(Pi)

    return (P, S)
