# Calcul matrice de Householder

def mat_householder(a):
    n = a.shape[0]
   
    e1 = np.zeros(n)
    e1[0]=1
    
    u = a - np.linalg.norm(a) * e1
    u = u.reshape(n,1)
    
    u1 = np.matmul(u,u.T)
    u2 = np.matmul(u.T,u)
   
    H = np.eye(n) - 2 * (u1/u2)
    
    return H

# Factorisation QR

def QR_Householder(A):
  n = A.shape[0]
  R = np.copy(A)
  Q = np.eye(n)
 
  for i in range(n):  
      a = R[i:,i]
      H = np.eye(n)
    
      if len(a)> 1 :
          H[i:,i:] = mat_householder(a)
          R = np.matmul(H,R)
          Q= np.matmul(Q,H)
      
  return Q,R

print("Matrice Q:")
print(Q)

print("Matrice R:")
print(np.round(R))

print("Matrice originale A:")
print(A)

print("Verification A = QR")
print(Q@R)


# Factorisation QR + application systèmes linéaires

def reso_QR_Householder(A,b):
    n = A.shape[0]
    R = np.copy(A)
    Q = np.eye(n)

    # 1) Factorisation QR   
    for i in range(n):
        a = R[i:,i]
        H = np.eye(n)
        
        if len(a)> 1 :
            H[i:,i:] = mat_householder(a)
            R = np.matmul(H,R) 
            Q= np.matmul(Q, H)
            
  # 2) Résolution système linéaire  
  X = np.zeros(n)    
    Y = Q.T @ b
    #  Résolution Rx = Y (algorithme de remontée)
    for i in range(n-1,-1,-1):
        X[i] = (Y[i] - np.dot(R[i,(i+1):],X[i+1:]))/R[i,i]
    
  return X
