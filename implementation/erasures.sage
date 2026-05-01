import random

def FieldToBits(x):
    """
    computes binary representation of GF(q) element in big-endian order
    """
    FF = x.parent()
    if FF.is_prime_field(): 
        q = FF.order()
        alpha = ceil(log(q, 2))
        val = int(x)
        binary = format(val, f"0{alpha}b")
        return binary
    else: 
        # for extension field, represent x as a degree < d polynomial over base field
        # and concatenate binary representations of its coefficients
        F = FF.prime_subfield()
        d = FF.degree()
        c1 = x.list()
        c = reversed(c1 + [F(0)] * (max(0, d - len(c1)))) # list of coefficients of x as a polynomial over base field
        return ''.join(FieldToBits(k) for k in c)

def MatrixToBits(M):
    """ 
    computes binary representation of matrix M over GF(q) using FieldToBits
    """
    Mb = []
    rows, cols = M.nrows(), M.ncols()
    for i in range(rows):
        bits = []
        for j in range(cols):
            bits.append(FieldToBits(M[i, j]))
        Mb.append(bits)   
    return Mb

def ErasedMatrix(M, p): 
    """ 
    generates bit erasures in matrix M over GF(q)
    """ 
    Mb = MatrixToBits(M)
    Me = []
    for row in Mb:
        r = []
        for string in row:
            s = ''.join('⊥' if random.random() < p else b for b in string)
            r.append(s)
        Me.append(r)
    return Me

def Counter(M):
    """ 
    counts erasures in matrix M and finds column with minimum erasure rate 
    """ 
    total = 0
    erased = 0
    l = len(M[0]) # amount of columns
    counter = [0] * l
    for j in range(l): # columns
        for i in range(len(M)): # rows
            string = M[i][j]
            total += len(string)
            erasures = string.count('⊥')
            erased += erasures
            counter[j] += erasures
    erasures_count = min(counter)
    best_col = counter.index(erasures_count)
    return f"{erased}/{total}", best_col, erasures_count # erasure rate, best column, erasures in best column

def BitsToField(x, FF): 
    """ 
    converts bit string into GF(q) element
    """ 
    if FF.is_prime_field(): 
        val = int(x, 2)
        y = FF(val)
        return y
    else: 
        F = FF.prime_subfield()
        alpha = ceil(log(F.order(), 2))
        d = FF.degree()
        c = [F(int(x[i*alpha:(i+1)*alpha], 2)) for i in range(d)]
        c.reverse()  
        return FF(c)
