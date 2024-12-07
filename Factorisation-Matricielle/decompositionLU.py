import numpy as np

# Décomposition LU 

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

# Décomposition LU et application à la résolution de systèmes linéaires

def facto_LU(A,b):
    n = A.shape[0]
    A1 = np.copy(A) # On définit une copie de A 
    
    # 1) Décomposition LU :
    for k in range(n):
        for i in range(k+1,n):
            if A1[k,k] == 0:
                return "La décomposition LU ne peut pas être calculée (pivot nul)"
            else :
                alpha =  A1[i,k] / A1[k,k] 
                for j in range(k+1, n):
                    A1[i,j] = A1[i,j] - alpha * A1[k,j] 
                A1[i,k] = alpha
                
    U = np.triu(A1)
    L = np.tril(A1,-1) + np.eye(n)
    
    # 2) Résolution système linéaire :  
    
        # Résolution Ly = b (algorithme de descente)
    Y = np.zeros(n)
    for i in range(n):
        Y[i]= (b[i] - np.dot(L[i,:i],Y[:i]))/L[i,i]
        
        # Résolution Ux = y (algorithme de remontée)
    X = np.zeros(n)
    for i in range(n-1,-1,-1):
        X[i] = (Y[i] - np.dot(U[i,(i+1):],X[i+1:]))/U[i,i]
    return X
