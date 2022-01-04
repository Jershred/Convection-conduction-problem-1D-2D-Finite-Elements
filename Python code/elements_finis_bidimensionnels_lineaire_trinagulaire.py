# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 15:09:30 2021

@author: Jérémy
"""


import numpy as np
import sympy as sp
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d  # Fonction pour la 3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


def tables(Lx,Ly,Nx,Ny,e):
    """ Description : Creation a partir du nombre de point et de la longueur du domaine la table de coord. globale des noeuds et la table de connexion pour un maillage triangulaire

        Donnees :   * Lx - Float : Longueur selon la direction x
                    * Ly - Float : Longueur selon la direction y
                    * Nx - Int   : Nbre point de discretisation selon x
                    * Ny - Int   : Nbre point de discretisation selon y

        Resultats : * Noeud : Table coord. gloable Noeuds
                    * Tbc : Table de connexion
    """
    nx = Nx - 1 # Nbre element sur x
    ny = Ny - 1 # Nbre element sur y


    lx = np.linspace(0,Lx,Nx)
    ly = np.linspace(0,Ly,Ny)
    Noeud = np.zeros((Nx*Ny,2))
    if e=='triangle':
        Tbc = np.zeros((2*nx*ny,3),dtype='int')
    elif e == 'carre':
        Tbc = np.zeros((nx*ny,4),dtype='int')
    

    Ne = 0
    Nn = 0
    i=0
    j=0
    compteur = 0

    while j < Ny-1:     # On se deplace sur les points sur y
        i = 0
        while i< Nx:  # On se deplace sur les points sur x
            if e == 'triangle':
                if 0<i and (Ne+1)%2 == 0:
                    A1=j*Nx+i
                    xA1 = lx[i]
                    yA1 = ly[j]
                    A2=(j+1)*Nx+i
                    xA2 = lx[i]
                    yA2 = ly[j+1]
                    A3=(j+1)*Nx+i-1
                    xA3 = lx[i-1]
                    yA3 = ly[j+1]

                elif 0<=i<=Nx-2:
                    A1=j*Nx+i
                    xA1 = lx[i]
                    yA1 = ly[j]
                    A2=j*Nx+i+1
                    xA2 = lx[i+1]
                    yA2 = ly[j]
                    A3=(j+1)*Nx+i
                    xA3 = lx[i]
                    yA3 = ly[j+1]
                    i=i+1

                elif (i+1)%Nx == 0:
                    break;

            elif e == 'carre':
                if (i+1)%Nx == 0:
                    break;
                else:
                    A1=j*Nx+i
                    xA1 = lx[i]
                    yA1 = ly[j]
                    A2=j*Nx+i+1
                    xA2 = lx[i+1]
                    yA2 = ly[j]
                    A3=(j+1)*Nx+i+1
                    xA3 = lx[i+1]
                    yA3 = ly[j+1]
                    A4=(j+1)*Nx+i
                    xA4 = lx[i]
                    yA4 = ly[j+1]
                    i=i+1

            if e == 'triangle':
                Tbc[Ne,0]=int(A1)
                Tbc[Ne,1]=int(A2)
                Tbc[Ne,2]=int(A3)
            elif e == 'carre':
                Tbc[Ne,0]=int(A1)
                Tbc[Ne,1]=int(A2)
                Tbc[Ne,2]=int(A3)
                Tbc[Ne,3]=int(A4)

            Noeud[A1,0] = xA1
            Noeud[A1,1] = yA1
            Noeud[A2,0] = xA2
            Noeud[A2,1] = yA2
            Noeud[A3,0] = xA3
            Noeud[A3,1] = yA3

            if e == 'carre':
                Noeud[A4,0] = xA4
                Noeud[A4,1] = yA4
            Ne = Ne + 1 # Numero de element
        j=j+1
    return (Tbc,Noeud)


def MatK (i,A1, A2, A3, Noeud):
    #i=0 pour x, i=1 pour y
    K=np.zeros((ne,ne),dtype="float64");
    A=np.zeros((ne,ne),dtype="float64");

    
    K[0,0]=(Noeud[A2,i]-Noeud[A3,i])**2;
    K[0,1]=(Noeud[A2,i]-Noeud[A3,i])*(Noeud[A3,i]-Noeud[A1,i]);
    K[0,2]=(Noeud[A2,i]-Noeud[A3,i])*(Noeud[A1,i]-Noeud[A2,i]);
    K[1,0]=K[0,1];
    K[1,1]=(Noeud[A3,i]-Noeud[A1,i])**2;
    K[1,2]=(Noeud[A3,i]-Noeud[A1,i])*(Noeud[A1,i]-Noeud[A2,i]);
    K[2,0]=K[0,2];
    K[2,1]=K[1,2];
    K[2,2]=(Noeud[A1,i]-Noeud[A2,i])**2;
    
    
    A[0,0]=1.0;
    A[0,1]=Noeud[A1,0];
    A[0,2]=Noeud[A1,1];
    A[1,0]=1.0;
    A[1,1]=Noeud[A2,0];
    A[1,2]=Noeud[A2,1];
    A[2,0]=1.0;
    A[2,1]=Noeud[A3,0];
    A[2,2]=Noeud[A3,1];
    
    detA=np.linalg.det(A);
    
    return (K*(1/(2*detA)));
    
def H1 (x,y,A1, A2, A3, Noeud):
    A=np.zeros((ne,ne),dtype="float64");
    
    
    A[0,0]=1.0;
    A[0,1]=Noeud[A1,0];
    A[0,2]=Noeud[A1,1];
    A[1,0]=1.0;
    A[1,1]=Noeud[A2,0];
    A[1,2]=Noeud[A2,1];
    A[2,0]=1.0;
    A[2,1]=Noeud[A3,0];
    A[2,2]=Noeud[A3,1];
    
    detA=np.linalg.det(A);
    
    return ((1/(2*detA))*((Noeud[A2,0]*Noeud[A3,1]-Noeud[A3,0]*Noeud[A2,1])+(Noeud[A2,1]-Noeud[A3,1])*x+(Noeud[A3,0]-Noeud[A2,0])*y))
    
def H2 (x,y,A1, A2, A3, Noeud):
    A=np.zeros((ne,ne),dtype="float64");
    
    A[0,0]=1.0;
    A[0,1]=Noeud[A1,0];
    A[0,2]=Noeud[A1,1];
    A[1,0]=1.0;
    A[1,1]=Noeud[A2,0];
    A[1,2]=Noeud[A2,1];
    A[2,0]=1.0;
    A[2,1]=Noeud[A3,0];
    A[2,2]=Noeud[A3,1];
    
    detA=np.linalg.det(A);
    
    return ((1/(2*detA))*((Noeud[A3,0]*Noeud[A1,1]-Noeud[A1,0]*Noeud[A3,1])+(Noeud[A3,1]-Noeud[A1,1])*x+(Noeud[A1,0]-Noeud[A3,0])*y))  
    
            
def H3 (x,y,A1, A2, A3, Noeud):
    A=np.zeros((ne,ne),dtype="float64");
    
    
    A[0,0]=1.0;
    A[0,1]=Noeud[A1,0];
    A[0,2]=Noeud[A1,1];
    A[1,0]=1.0;
    A[1,1]=Noeud[A2,0];
    A[1,2]=Noeud[A2,1];
    A[2,0]=1.0;
    A[2,1]=Noeud[A3,0];
    A[2,2]=Noeud[A3,1];
    
    detA=np.linalg.det(A);
    
    return ((1/(2*detA))*((Noeud[A1,0]*Noeud[A2,1]-Noeud[A2,0]*Noeud[A1,1])+(Noeud[A1,1]-Noeud[A2,1])*x+(Noeud[A2,0]-Noeud[A1,0])*y))

def MatMx (x,y,A1, A2, A3, Noeud,b1,b2):
    M=np.zeros((ne,ne),dtype="float64");
    
    M[0,0]=sp.integrate(H1(x,y,A1, A2, A3, Noeud)**2,(x,b1,b2));
    M[0,1]=sp.integrate(H1(x,y,A1, A2, A3, Noeud)*H2(x,y,A1, A2, A3, Noeud),(x,b1,b2));
    M[0,2]=sp.integrate(H1(x,y,A1, A2, A3, Noeud)*H3(x,y,A1, A2, A3, Noeud),(x,b1,b2));
    M[1,0]=M[0,1]
    M[1,1]=sp.integrate(H2(x,y,A1, A2, A3, Noeud)**2,(x,b1,b2));
    M[1,2]=sp.integrate(H3(x,y,A1, A2, A3, Noeud)*H2(x,y,A1, A2, A3, Noeud),(x,b1,b2));
    M[2,0]=M[0,2]
    M[2,1]=M[1,2]
    M[2,2]=sp.integrate(H3(x,y,A1, A2, A3, Noeud)**2,(x,b1,b2));
    
  
    return (M)

def MatMy (x,y,A1, A2, A3, Noeud,b1,b2):
    M=np.zeros((ne,ne),dtype="float64");
    
    M[0,0]=sp.integrate(H1(x,y,A1, A2, A3, Noeud)**2,(y,b1,b2));
    M[0,1]=sp.integrate(H1(x,y,A1, A2, A3, Noeud)*H2(x,y,A1, A2, A3, Noeud),(y,b1,b2));
    M[0,2]=sp.integrate(H1(x,y,A1, A2, A3, Noeud)*H3(x,y,A1, A2, A3, Noeud),(y,b1,b2));
    M[1,0]=M[0,1]
    M[1,1]=sp.integrate(H2(x,y,A1, A2, A3, Noeud)**2,(y,b1,b2));
    M[1,2]=sp.integrate(H3(x,y,A1, A2, A3, Noeud)*H2(x,y,A1, A2, A3, Noeud),(y,b1,b2));
    M[2,0]=M[0,2]
    M[2,1]=M[1,2]
    M[2,2]=sp.integrate(H3(x,y,A1, A2, A3, Noeud)**2,(y,b1,b2));
    
  
    return (M)

def VectBy(x,y,A1, A2, A3, Noeud,b1,b2):
    B=np.zeros((ne,1),dtype="float64");
    
    B[0]=sp.integrate(H1(x,y,A1, A2, A3, Noeud),(y,b1,b2));
    B[1]=sp.integrate(H2(x,y,A1, A2, A3, Noeud),(y,b1,b2));
    B[2]=sp.integrate(H3(x,y,A1, A2, A3, Noeud),(y,b1,b2));
    
    return (B*(h/lbda)*Tinf)

def VectBx (x,y,A1, A2, A3, Noeud,b1,b2):
    B=np.zeros((ne,1),dtype="float64");
    
    B[0]=sp.integrate(H1(x,y,A1, A2, A3, Noeud),(x,b1,b2));
    B[1]=sp.integrate(H2(x,y,A1, A2, A3, Noeud),(x,b1,b2));
    B[2]=sp.integrate(H3(x,y,A1, A2, A3, Noeud),(x,b1,b2));
    
    return (B*(h/lbda)*Tinf)



#MAIN


#-----Paramétres------

H=4.0;  #en m  
r=2.0;  #en m

To=100.0;   #en °C
Tinf=20.0;  #en °C
h=80.0;     #en W/m*2*°C
lbda=40.0;  #en W/m*°C


NE=15;  #nombre de noeuds
ne=3;   #taille matrice de base
Nx=3;   #Nombre point de discretisation selon x
Ny=5;  #Nombre point de discretisation selon y

#----Creation maillage---
(Tbc,Noeud)=tables (r,H,Nx,Ny,'triangle')


x, y = sp.symbols('x y');

#-----Assemblage-------


A=np.zeros((NE,NE),dtype="float64");    #matrice 
B=np.zeros((NE,1),dtype="float64");     #second membre
Kx=np.zeros((ne,ne),dtype="float64");
Ky=np.zeros((ne,ne),dtype="float64");
M1=np.zeros((ne,ne),dtype="float64");
b=np.zeros((ne,1),dtype="float64");




#Assemblage avec les matrices Kx et Ky:
    
for k in range(NE+1):
    list=np.array([Tbc[k,0],Tbc[k,1],Tbc[k,2]]);
    
    Kx=MatK (0,list[0], list[1], list[2], Noeud);
    Ky=MatK (1,list[0], list[1], list[2], Noeud);
    
    # print('Kx:\n',Kx)
    # print('Ky:\n',Ky)
    
   
    
    for i in range(ne):
        for j in range(ne):
            A[list[i],list[j]]= A[list[i],list[j]]+Kx[i,j]+Ky[i,j];
            
            
            
           
#Assemblage avec M et B pour les éléments sur bords impliquant une condidtion de Neumann:
 
    
#Sur le bord à x=r:
    
E_list=np.array([3,7,11,15]);   #liste des éléments sur le bord
ylist=np.linspace(0,H,5);   #intervalle d'intégration

for k in range(np.size(ylist)-1):
    
    list=np.array([Tbc[E_list[k],0], Tbc[E_list[k],1], Tbc[E_list[k],2]]);#table avec les noeuds d'un élément k
    
    
    M1=MatMy(r, y, list[0], list[1], list[2], Noeud, ylist[k], ylist[k+1]);
    
    b=VectBy(r, y, list[0], list[1], list[2], Noeud, ylist[k], ylist[k+1]);
    # print('b1 et b2:\n',ylist[k],ylist[k+1])
    print('My élémentaire:\n',M1)
    print('by élémentaire:\n',b)
    
             
    for i in range(ne):
        B[list[i]]=B[list[i]]+b[i];
       
        for j in range(ne):
            A[list[i],list[j]]= A[list[i],list[j]]+M1[i,j];  
 

#Sur le bord à y=H:  
    
E_list1=[13,15];#liste d'éléments sur le bord
xlist=np.linspace(0,r,3);#intervalle d'intégration

for k in range(np.size(xlist)-1):
    
    list=np.array([Tbc[E_list1[k],0], Tbc[E_list1[k],1], Tbc[E_list1[k],2]]);#table avec les noeud d'un élément k
    
    M1=MatMx(x, H,  list[0], list[1], list[2], Noeud, xlist[k], xlist[k+1]);
   
    b=VectBx(x, H, list[0], list[1], list[2], Noeud, xlist[k], xlist[k+1]);
    
    print('Mx élémentaire:\n',M1)
    print('bx élémentaire:\n',b)
        
    for i in range(ne):
        B[list[i]]=B[list[i]]+b[i];
        for j in range(ne):
            A[list[i],list[j]]= A[list[i],list[j]]+M1[i,j];   


#----Conditions de dirichlet----

E_bord=np.array([0,1,2]);#elements sur les bords où T(0)=To

for i in range(NE):
    for p in range(np.size(E_bord)):
      if i==E_bord[p]: 
        A[i,:]=0.0;
        A[i,i]=1.0
        B[i,0]=To;

#---résolution systéme----

T=np.linalg.solve(A,B);


#Tracé


# Tableau pour les 3 axes

x=np.linspace(0,r,Nx);
y=np.linspace(0,H,Ny);
X, Y= np.meshgrid(x, y)

# Z=np.zeros(NE)
# for i in range(NE):
#     Z[i]=T[i,0]

Z=np.array([[T[0,0],T[1,0],T[2,0]],
            [T[3,0],T[4,0],T[5,0]],
            [T[6,0],T[7,0],T[8,0]],
            [T[9,0],T[10,0],T[11,0]],
            [T[12,0],T[13,0],T[14,0]]])

print('Z;\n',Z)


figure(1)

gca(projection='3d').plot_surface(X,Y,Z,cmap=cm.coolwarm, linewidth=0)
xlabel('X(m)')
ylabel('Y(m)')
title('Température dans la demi aillette')

#Courbes de niveau:
figure(2)

cont=contour(X,Y,Z,cmap=cm.nipy_spectral_r)
clabel(cont,fmt='%d')
xlabel('X (m)')
ylabel('Y (m)')
title('Tracé de contour de la température')