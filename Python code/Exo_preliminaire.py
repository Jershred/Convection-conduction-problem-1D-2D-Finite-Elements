# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 12:39:11 2021

@author: Jershred
"""

import numpy as np
import matplotlib.pyplot as plt
from random import*

def pentadiagonale(n):
    M=np.zeros((n**2,n**2))
    for i in range (n**2):
        for j in range (n**2):
            if i==j:
                M[i,j]=4
            if i==j+1:
                M[i,j]=-1
            if i==j-1:
                M[i,j]=-1
            if i==j-n:
                M[i,j]=2
            if i==j+n:
                M[i,j]=2
    return(M)
    
def sousblock(n,A):
    """A=np.array([[2,1,1],[1,2,1],[1,1,2]])"""
    a=np.zeros((9,9))
    for k in range(0,7,n):
        for i in range(3):
            for j in range(3):
                a[i+k,j+k]=a[i+k,j+k]+A[i,j]
    return(a)

"""x = np.array([randrange(100),randrange(100),randrange(100),randrange(100)])
y = np.array([randrange(100),randrange(100),randrange(100),randrange(100)])
plt.scatter(x, y, label="point random")
plt.xlabel("abscisses")
plt.ylabel("ordonnees")
plt.yscale('log')
plt.grid(True,which="both", linestyle='--')
plt.legend()
plt.show()"""

def factorisation(A):
    n=len(A)
    U=np.zeros(n,n)
    L=np.zeros(n,n)
    for v in range(n):
        for w in range(n):
            if v==w:
                L[v,w]=1
    for j in range(n):
        for i in range(j):
            S1=0
            for k in range(i):
                S1=S1+L[i,k]*U[k,j]
            U[k,l]=A[k,l]-S1
        for i in range(j,n):
            S2=0
            for k in range(j):
                S1=S1+L[i,k]*U[k,j]
            L[i,j]=1/U[j,j]*(A[i,j]-S2)
            
def tables(Lx=10,Ly=20,Nx=3,Ny=5,e='triangle'):
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

# if __name__ == "__main__":
#     [Tbc,Coord]= tables(e='carre')
#     print (Tbc)
#     print (np.shape(Tbc))
#     print (Coord)
#     maillage.maillage(Tbc,Coord)

#     plt.show()

def triangle(x1,y1,x2,y2,x3,y3):
    plt.plot([x1, x2], [y1, y2], 'r-', lw=2) # Red straight line
    plt.plot([x2, x3], [y2,y3], 'r-', lw=2)
    plt.plot([x3, x1], [y3,y1], 'r-', lw=2)
    plt.show()
    return()

def trace_maillage(Lx=10,Ly=20,Nx=3,Ny=5,e='triangle'):
    [Tbc,Coord]= tables(e='carre')
    for i in range(len(Tbc)):
        triangle(Coord(i,0),Coord(i,1),Coord(i+1,0),Coord(i+1,1),Coord(i+2,0),Coord(i+2,1))
    return()
    
    
    
    
    
    
    
    
    