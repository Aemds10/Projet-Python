import numpy as np


def facto_LU(A):
    n = A.shape[0]
    A1 = np.copy(A)
    
    for k in range(n):
        for i in range(k+1,n):
            if A1[k,k] == 0:
                return "La décomposition LU ne peut pas être calculée (pivot nul)"
            else :
                alpha =  A1[i,k] / A1[k,k] # déf alpha
                for j in range(k+1, n):
                    A1[i,j] = A1[i,j] - alpha * A1[k,j]  # a_ij => U
                 A1[i,k] = alpha  # => L  
    return A1

# Matrice U
U = np.triu(resA)
print("Matrice U:")
print(U)

# Matrice L
L = np.tril(resA,-1) + np.eye(resA.shape[0])
print("Matrice L:")
print(L)
