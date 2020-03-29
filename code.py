import numpy as np

eps=0.5e-5
A1=np.array([[2.,4.,6.],[-2.,1.,1.],[-1.,-1.,2.]])
A2=np.array([[1.,1+0.5e-15,3.],[2.,2.,20.],[3.,6.,4.]])
A3=np.array([[1.,1+eps,3.],[2.,2.,20.],[3.,6.,4.]])
X=np.zeros((3,1))
B1=np.array([6.,-1.,-2.])
B2=np.array([5.+0.5e-15,24.,13.])
B3=np.array([5+eps,24.,13.])



print("*********************MATRICE A1*************************")
#factorisation LU de la matrice A1
U=np.zeros((3,3))
#initialisation de la matrice triangulaire superieure
L=np.eye(3)
#initialisation de la matrice triangulaire inferieure a diagonale unite
for k in range(3):
	for c in range(k,3):
		U[k,c]=A1[k,c]
		for j in range(k):
			U[k,c]=U[k,c]-L[k,j]*U[j,c]
	for d in range(k+1,3):
		L[d,k]=A1[d,k]
		for j in range(k):
			L[d,k]=L[d,k]-L[d,j]*U[j,k]
		L[d,k]=L[d,k]/U[k,k]
print("factorisation LU de A1")
print("L =")
print(L)
print("U =")
print(U)

#resolution du systeme
#methode 1:
#resolution de Ly=b
y=np.zeros((3,1))
for i in range(3):
	y[i]=B1[i]
	for j in range (i):
		y[i]=y[i]-L[i,j]*y[j]
	y[i]=y[i]/L[i,i]
#resolution de Ux=y
for i in range(2,-1,-1):
	X[i]=y[i]	
	for j in range(i+1,3):
		X[i]=X[i]-U[i,j]*X[j]
	X[i]=X[i]/U[i,i]
print("La solution de methode 1 est : ")
print(X)

#methode 2:
A0=np.linalg.inv(L)@L@U
B0=np.linalg.inv(L)@B1
for i in range(2,-1,-1):
	X[i]=B0[i]
	for j in range(i+1,3):
		X[i]=X[i]-X[j]*A0[i,j]
	X[i]=X[i]/A0[i,i]
print("La solution de methode 2 est : ")
print(X)


#factorisation avec pivotage de matrice A1
#méthode 1
A0=A1
E1=np.eye(3)
for i in range(2):
	E=np.eye(3)
	for j in range(i+1,3):
		E[j,i]=-A0[j,i]/A0[i,i]
	E1=E@E1
	A0=E@A0
print("factorisation avec pivotage de A1")
print("L =")
print(L)
print("U =")
print(A0)

#méthode 2
A0=A1
L=np.eye(3)
for i in range(2):
	for j in range(i+1,3):
		L[j,i]=A0[j,i]/A0[i,i]
		for k in range(i,3):
			A0[j,k]=A0[j,k]-L[j,i]*A0[i,k]
print("factorisation avec pivotage de A1")
print("L =")
print(L)
print("U =")
print(A0)

#resolution du systeme
B0=np.linalg.inv(L)@B1
for i in range(2,-1,-1):
	X[i]=B0[i]
	for j in range(i+1,3):
		X[i]=X[i]-X[j]*A0[i,j]
	X[i]=X[i]/A0[i,i]
print("La solution est :")
print(X)
print(' ')



print("*********************MATRICE A2*************************")
#factorisation LU de matrice A2
U=np.zeros((3,3))
#initialisation de la matrice triangulaire superieure
L=np.eye(3)
#initialisation da la matrice triangulaire inferieure a diagonale unite
for k in range(3):
	for c in range(k,3):
		U[k,c]=A2[k,c]
		for j in range(k):
			U[k,c]=U[k,c]-L[k,j]*U[j,c]
	for d in range(k+1,3):
		L[d,k]=A2[d,k]
		for j in range(k):
			L[d,k]=L[d,k]-L[d,j]*U[j,k]
		L[d,k]=L[d,k]/U[k,k]
print("factorisation LU de A2")
print("L =")
print(L)
print("U =")
print(U)

#resolution du systeme
#resolution de Ly=b
y=np.zeros((3,1))
for i in range(3):
	y[i]=B2[i]
	for j in range (i):
		y[i]=y[i]-L[i,j]*y[j]
	y[i]=y[i]/L[i,i]
#resolution de Ux=y
for i in range(2,-1,-1):
	X[i]=y[i]	
	for j in range(i+1,3):
		X[i]=X[i]-U[i,j]*X[j]
	X[i]=X[i]/U[i,i]
print("La solution de methode 1 est : ")
print(X)


#factorisation avec pivotage de matrice A2
A0=A2
E1=np.eye(3)
for i in range(2):
	E=np.eye(3)
	for j in range(i+1,3):
		E[j,i]=-A0[j,i]/A0[i,i]
	E1=np.dot(E,E1)
	A0=np.dot(E,A0)
	#A0=E(3)*E(2)*A2(1)=E1*A2(1)
print("factorisation avec pivotage de A2")
print("L =")
print(E1)
print("U =")
print(A0)

#resolution du systeme
B0=E1@B2
for i in range(2,-1,-1):
	X[i]=B0[i]
	for j in range(i+1,3):
		X[i]=X[i]-X[j]*A0[i,j]
	X[i]=X[i]/A0[i,i]
print("La solution est : ")
print(X)
print(' ')



print("*********************MATRICE A3*************************")
#factorisation LU de matrice A3
U=np.zeros((3,3))
#initialisation de la matrice triangulaire superieure
L=np.eye(3)
#initialisation da la matrice triangulaire inferieure a diagonale unite
for k in range(3):
	for c in range(k,3):
		U[k,c]=A3[k,c]
		for j in range(k):
			U[k,c]=U[k,c]-L[k,j]*U[j,c]
	for d in range(k+1,3):
		L[d,k]=A3[d,k]
		for j in range(k):
			L[d,k]=L[d,k]-L[d,j]*U[j,k]
		L[d,k]=L[d,k]/U[k,k]
print("factorisation LU de A3")
print("L =")
print(L)
print("U =")
print(U)

#resolution du systeme
#resolution de Ly=b
y=np.zeros((3,1))
for i in range(3):
	y[i]=B3[i]
	for j in range (i):
		y[i]=y[i]-L[i,j]*y[j]
	y[i]=y[i]/L[i,i]
#resolution de Ux=y
for i in range(2,-1,-1):
	X[i]=y[i]	
	for j in range(i+1,3):
		X[i]=X[i]-U[i,j]*X[j]
	X[i]=X[i]/U[i,i]
print("La solution de methode 1 est : ")
print(X)


#factorisation avec pivotage de matrice A3
A0=A3
E1=np.eye(3)
for i in range(2):
	E=np.eye(3)
	for j in range(i+1,3):
		E[j,i]=-A0[j,i]/A0[i,i]
	E1=np.dot(E,E1)
	A0=np.dot(E,A0)
	#A0=E(3)*E(2)*A3(1)=E1*A3(1)
print("factorisation avec pivotage de A3")
print("L =")
print(E1)
print("U =")
print(A0)

#resolution du systeme
B0=E1@B3
for i in range(2,-1,-1):
	X[i]=B0[i]
	for j in range(i+1,3):
		X[i]=X[i]-X[j]*A0[i,j]
	X[i]=X[i]/A0[i,i]
print("La solution est : ")
print(X)
