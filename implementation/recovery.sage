"""
for one_vector_to_key and in_secret_subspace (these methods were introduced by Pierre Pébereau in https://doi.org/10.1007/978-3-031-62746-0_5), 
we used code https://github.com/pi-r2/OneVector, which is available under MIT license:

The MIT License (MIT)

Copyright (c) 2026 Pierre Pébereau

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

def one_vector_to_key(P, v, o):
    """ 
    This function takes as input a UOV public key P and vector v in Oil subspace and returns basis of O 
    """
    charac = P[0].base_ring().characteristic()
    if charac == 2 :
        J = matrix([v*(p + p.transpose()) for p in P])
    else:
        J = matrix([v*p for p in P])
    B = matrix(J.right_kernel().basis())
    B2 = []
    for p in P :
        phat = B*p*B.transpose() #Restriction to ker(J(v))
       	if charac == 2 :
            phat = phat + phat.transpose()
        for b in phat.kernel().basis() :
            if len(B2) == 0 or b not in span(B2) :
                B2.append(b)
        if len(B2) == o :
            break
    B3 = matrix(B2)
        
    C = B3*B
    return C

def in_secret_subspace(P, v, o):
    """
    This function takes as input a UOV public key P, vector v and Oil subspace dimension o 
    and returns True if the vector belongs to the secret subspace, and False otherwise. 
    """
    charac = P[0].base_ring().characteristic()

    for p in P : #Sanity check
        if v*p*v != 0 :
            return False
    
    if charac == 2 :
        J = matrix([v*(p + p.transpose()) for p in P])
    else:
        J = matrix([v*p for p in P])

    B = matrix(J.right_kernel().basis()) # B is the basis of ker(J(v))
    d = J.right_kernel().dimension()

    if d < o:
        return False
    
    for p in P :
        phat = B*p*B.transpose() #Restriction to K(v)
        if charac == 2 :
            phat = phat + phat.transpose()
        if phat.rank() > 2*(d-o) : 
            return False 
    return True 
