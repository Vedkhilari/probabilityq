import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.image as mpimg
import math
sys.path.insert(0, '/home/void/randomvector/codes/CoordGeo/')
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen
#If using termux
import subprocess
import shlex
#end if

#Sample size
simlen = 2
#Possible outcomes
#n = range(2,13)
# Generate X1 and X2
#y = np.random.randint(-6,6, size=(3, simlen))
y = np.random.randint(-6,6, size=(3, simlen))
print(y)
A=y[0]
B=y[1]
C=y[2]

#1.1.1

d = B- A
e = C - B
f = A - C

print("The direction vector of AB is ",d)
print("The direction vector of BC is ",e)
print("The direction vector of CA is ",f)

#1.1.2
BC_matrix = B-C

length_BC = np.linalg.norm(BC_matrix)

print("Length of side BC:", length_BC)

AB_matrix = B-A 

length_AB = np.linalg.norm(AB_matrix)

print("Length of side AB:", length_AB)

CA_matrix = C-A

length_CA = np.linalg.norm(CA_matrix)

print("Length of side CA:", length_CA)

#1.1.3
Mat = np.array([[1,1,1],[A[0],B[0],C[0]],[A[1],B[1],C[1]]])

rank = np.linalg.matrix_rank(Mat)

if (rank<=2):
	print("Hence proved that points A,B,C in a triangle are collinear")
else:
	print("The given points are not collinear")
	
	
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
#plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
#plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
#D = D.reshape(-1,1)
#E = E.reshape(-1,1)
#F = F.reshape(-1,1)
#G = G.reshape(-1,1)
tri_coords = np.block([[A, B, C]])
plt.scatter(tri_coords[0, :], tri_coords[1, :])
vert_labels = ['A', 'B', 'C']
for i, txt in enumerate(vert_labels):
    offset = 10 if txt == 'C' else -10
    plt.annotate(txt,
                 (tri_coords[0, i], tri_coords[1, i]),
                 textcoords="offset points",
                 xytext=(0, offset),
                 ha='center')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

#if using termux
#plt.savefig('tri_sss.pdf')
plt.savefig('/home/void/randomvector/figure1.png')


#1.1.4
m1=(B-A)
m2=(C-B)
m3=(A-C)
print("parametric of AB form is x:",A,"+ k",m1)
print("parametric of BC form is x:",B,"+ k",m2)
print("parametric of CA form is x:",C,"+ k",m3)

#1.1.5

#getting the equation of line
omat = np.array([[0,1],[-1,0]])
m = B-A   #direction vector
n = omat@m    #normal vector
c = n@A
eqn = f"{n}x = {c}"
print("The equation of line AB is",eqn)

#plotting the line AB
x_AB = line_gen(A, B)
A = A.reshape(-1,1)
B = B.reshape(-1,1)
tri_coords = np.block([A,B])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
vert_labels=['A','B']     #for labelling points A and B 
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid()
plt.savefig('/home/void/randomvector/figure2.png')


#1.1.6
def AreaCalc(A, B, C):
    AB = A - B
    AC = A - C
#cross_product calculation
    cross_product = np.cross(AB,AC)
#magnitude calculation
    magnitude = np.linalg.norm(cross_product)

    area = 0.5 * magnitude

    return area
area_ABC = AreaCalc(A, B, C)
print("Area of triangle ABC:", area_ABC)

#1.1.7

dotA=((B-A).T)@(C-A)
dotA=dotA[0,0]
NormA=(np.linalg.norm(B-A))*(np.linalg.norm(C-A))
print('value of angle A: ', np.degrees(np.arccos((dotA)/NormA)))


dotB=(A-B).T@(C-B)
dotB=dotB[0,0]
NormB=(np.linalg.norm(A-B))*(np.linalg.norm(C-B))
print('value of angle B: ', np.degrees(np.arccos((dotB)/NormB)))

dotC=(A-C).T@(B-C)
dotC=dotC[0,0]
NormC=(np.linalg.norm(A-C))*(np.linalg.norm(B-C))
print('value of angle C: ', np.degrees(np.arccos((dotC)/NormC)))



#1.2.1
D = (B + C)/2

#Similarly for E and F
E = (A + C)/2
F = (A + B)/2

print("D:", list(D))
print("E:", list(E))
print("F:", list(F))

#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)


#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')


#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
tri_coords = np.block([[A,B,C,D,E,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','E','F']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

#if using termux
#plt.savefig('tri_sss.pdf')
plt.savefig('/home/void/randomvector/figure3.png')

#1.2.2
D = (B + C) / 2
E = (C + A) / 2
F = (A + B) / 2

def dir_vec(A, B):
    return B - A

def norm_vec(A, B):
    return omat @ dir_vec(A, B) 

def line_gen(A, B):
    len = 10
    dim = A.shape[0]
    x_AB = np.zeros((dim, len))
    lam_1 = np.linspace(0, 1, len)
    for i in range(len):
        temp1 = A + lam_1[i] * (B - A)
        x_AB[:, i] = temp1.T
    return x_AB


x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)

plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')

x_AD = line_gen(A, D)
plt.plot(x_AD[0, :], x_AD[1, :], label='$AD$')

x_BE = line_gen(B, E)
plt.plot(x_BE[0, :], x_BE[1, :], label='$BE$')

x_CF = line_gen(C, F)
plt.plot(x_CF[0, :], x_CF[1, :], label='$CF$')

A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
tri_coords = np.block([[A,B,C,D,E,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','E','F']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-10,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')


plt.show()



#1.2.3
mat = np.array([[0,1],[-1,0]]) 

np.set_printoptions(precision=2)

def midpoint(x,y):
	return (x+y)/2

def norm_vec(A,B):
	return np.matmul(omat, dir_vec(A,B))

#direction
def dir_vec(A,B):
	return B-A
	
#Given D,E,F are midpoints of BC,CA,AB
D=midpoint(B,C)
E=midpoint(C,A)
F=midpoint(A,B)

#intersection of lines
def line_intersect(n1,A1,n2,A2):
	N=np.block([[n1],[n2]])
	p = np.zeros(2)
	p[0] = n1@A1
	p[1] = n2@A2
	#Intersection
	P=np.linalg.inv(N)@p
	return P
G=line_intersect(norm_vec(F,C),C,norm_vec(E,B),B)
print("("+str(G[0])+","+str(G[1])+")")

#Hence verified that A - F = E - D and AFDE is a parallelogram

#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_BE = line_gen(B,E)
x_CF = line_gen(C,F)


#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
G = G.reshape(-1,1)
tri_coords = np.block([[A, B, C, D, E, F, G]])
plt.scatter(tri_coords[0, :], tri_coords[1, :])
vert_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
for i, txt in enumerate(vert_labels):
    offset = 10 if txt == 'G' else -10
    plt.annotate(txt,
                 (tri_coords[0, i], tri_coords[1, i]),
                 textcoords="offset points",
                 xytext=(0, offset),
                 ha='center')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

#if using termux
#plt.savefig('tri_sss.pdf')
plt.savefig('/home/void/randomvector/figure4.png')


#1.2.4
#initializing variables
G = np.zeros(2) 
D = np.zeros(2)
E = np.zeros(2)
F = np.zeros(2)

#find D,E,F
D = line_section(B,C,1.0)
E = line_section(C,A,1.0)
F = line_section(A,B,1.0)

G = line_intersect(B,E,C,F)

print("The vector G is",G)

#finding the direction of vector
BG = dir_vec(B,G)
GE = dir_vec(G,E)
GF = dir_vec(G,F)
CG = dir_vec(C,G)
AG = dir_vec(A,G)
GD = dir_vec(G,D)

#finding the norm of vector
n_BG = norm_vec(B,G)
n_GE = norm_vec(G,E)
n_GF = norm_vec(G,F)
n_CG = norm_vec(C,G)
n_AG = norm_vec(A,G)
n_GD = norm_vec(G,D)

print("The direction vector of BG is",BG)
print("The direction vector of GE is",GE)
print("The direction vector of GF is",GF)
print("The direction vector of CG is",CG)
print("The direction vector of AG is",AG)
print("The direction vector of GD is",GD)

print("The norm of the vector BG is",n_BG)
print("The norm of the vector GE is",n_GE)
print("The norm of the vector GF is",n_GF)
print("The norm of the vector CG is",n_CG)
print("The norm of the vector AG is",n_AG)
print("The norm of the vector GD is",n_GD)

print("The ratio BG/GE is",n_BG/n_GE)
print("The ratio CG/GF is",n_CG/n_GF)
print("The ratio AG/GD is",n_AG/n_GD)

#1.2.5
D=(B+C)/2

print("The mid point of B and C is {D}")

G=(A+B+C)/3

print("The centroid of triangleABC is {G}")

Mat = np.array([[1,1,1],[A[0],D[0],G[0]],[A[1],D[1],G[1]]])

rank = np.linalg.matrix_rank(Mat)

if (rank==2):
	print("Hence proved that points A,G,D in a triangle are collinear")
else:
	print("Error")

#1.2.6
G = (A + B + C) / 3
print("centroid of the given triangle: ")      
      
print(G)
     
print("Hence Q.1.2.6 is verified.")

#1.2.7
#Now we calculate the co-ordinates of the mid-points D,E,F of the sides AB,BC,CA respectively

D = (B + C)/2
E = (A + C)/2
F = (A + B)/2

print(f"A - F = {A-F}")
print(f"E - D = {E-D}")

#Hence verified that A - F = E - D and AFDE is a parallelogram

#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_DE = line_gen(D,E)
x_DF = line_gen(D,F)


#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_DE[0,:],x_DE[1,:],label='$DE$')
plt.plot(x_DF[0,:],x_DF[1,:],label='$DF$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
tri_coords = np.block([[A,B,C,D,E,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','E','F']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

#if using termux
#plt.savefig('tri_sss.pdf')
plt.savefig('/home/void/randomvector/figure5.png')





