'''
here, we use msolve library https://msolve.lip6.fr to solve quadratic systems over prime fields
and native sage to solve quadratic systems over extension fields
'''
load("uov.sage")
load("recovery.sage")
load("erasures.sage")
import itertools
set_verbose(-1)

def sim(v, o, m, q, p, wtilde):
    '''
    Simulation of partial key exposure attack in bit erasure model on UOV-schemes
    p is bit-erasure probability, wtilde is enumeration parameter
    '''
    print('v = ', v, ' o = ', o, ' m = ', m, 'q = ', q)
    FF = GF(q)
    n = v + o
    print('Key Generation...')
    P, S = KeyGen(v, o, m, q)

    Sinv = S.inverse()
    print("\nOil subspace basis in column order (n × o):")
    for row in Sinv[:n, v:]:
        print(list(row))
    
    B = Sinv[:v, v:]
    print()
    print("matrix B (v × o):")
    for i in range(v):
        print(list(B.row(i)))

    SinvE = ErasedMatrix(B, p)
    print("\nErased B:")
    for row in SinvE:
        print("[" + " ".join(row) + "]")

    rate, BestColumnIndex, erasures_best_col = Counter(SinvE) # erasure rate, best column index, erasures in best column
    BestColumn = [row[BestColumnIndex] for row in SinvE]
    BestColumnString = " ".join(BestColumn)

    variables = [[BestColumn[i].count('⊥'), i, BestColumn[i]] for i in range(len(BestColumn))]
    variables = [var for var in variables if var[0] > 0] # pick only erased variables
    variables.sort() # sort them

    w = len(variables)

    if w - wtilde < 0 or w - wtilde > m:
        print(f"\nAttack failed:" f"\nexpected wtilde ≤ w = {w} and w - wtilde ≤ m = {m}"
             f"\ngot wtilde = {wtilde}, w - wtilde = {w-wtilde}")
        return

    print(f"\nBest column (index = {BestColumnIndex}, erased coordinates = {w}, erasures = {erasures_best_col}):")
    print(BestColumnString)

    Wtilde = variables[:wtilde] # pick first "wtilde" of sorted variables

    print(f"We took wtilde = {wtilde}, Wtilde variables:")
    show = [f"{string} (index = {index}, erasures = {count})" for count, index, string in Wtilde]
    print(" ".join(show))
    print("\nRunning attack")

    R = PolynomialRing(FF, n, 'x', order='degrevlex')
    x = vector(R.gens())
    eq = [(x * pi * x) for pi in P] # x^T*P*x = 0

    solution = None
    success = False
    enum = sum(count for count, _, _ in Wtilde) # amount of bits for enumeration step

    ErasurePositions = []
    for count, index, string in Wtilde:
        positions = [position for position, element in enumerate(string) if element == '⊥']
        ErasurePositions.append((index, positions)) # variable's index, positions of erasures

    e_j = [0] * o
    e_j[BestColumnIndex] = 1
    print(f"Guessing {enum} bits")
    print(f"This might take some time...")

    for bits in itertools.product([0, 1], repeat=enum):
        constraints = []    
        for i in range(o):
            constraints.append(x[v + i] - FF(e_j[i])) # we look for j-th Oil subspace basis vector

        counter = 0 # bit counter
        for i in range(v):
            string = BestColumn[i]       
            if '⊥' not in string: # variable has no erasures
                constraints.append(x[i] - BitsToField(string, FF)) 
            else:
                InWtilde = False
                for idx, positions in ErasurePositions:
                    if i == idx: # if variable is in Wtilde
                        temp = list(string)
                        for pos in positions:
                            temp[pos] = str(bits[counter])
                            counter += 1
                        constraints.append(x[i] - BitsToField("".join(temp), FF))
                        InWtilde = True
                        break           
                if not InWtilde:
                    continue

        I = R.ideal(eq + constraints)
        try:
            sols = I.variety(algorithm='msolve', proof=False) if FF.is_prime_field() else I.variety()
            if sols:
                for s in sols:
                    res = vector(FF, [s[var] for var in x])
                    bool1 = all((res * pi * res == 0) for pi in P) # check if res lies in forgery variety
                    bool2 =  in_secret_subspace(P, res, o) # check if res lies in Oil subspace
                    print(f"{res}\nP(solution) = 0: {bool1}, solution in Oil subspace: {bool2}")
                    if in_secret_subspace(P, res, o):
                        print(f"Obtained oil vector: {res}")
                        basis = one_vector_to_key(P, res, o)
                        if basis.row_space() == span(S.inverse().columns()[v:]):
                            solution = res
                            success = True
                            print("Oil subspace recovered")
                            break
                        else:
                            print("Failed to recover Oil subspace")   
        except Exception:
            continue
        if success:
            break
    if not success:
        print("\nAttack failed")

print("=" * 150)
print("\nTest 1 (toy QR-UOV attack)") 
sim(20, 14, 31, 7, 0.8, 6) # v, o, m, q, p, wtilde
print("=" * 150)
print("\nTest 2 (toy MAYO attack)")
sim(15, 10, 25, 8, 0.8, 6)
print("=" * 150)
print("\nTest 3 (toy UOV attack)")
sim(15, 10, 10, 7, 0.8, 6) 

# In test 3, solution space contains approximately 7^(15) vectors, and we have 7^(10) additional zeros in Oil subspace,
# that's why we may obtain a lot of false vectors, which are solutions to P(x) = 0, but don't lie in Oil subspace
