# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 15:21:14 2021

@author: Jershred
"""

import numpy as np
import matplotlib.pyplot as plt
from random import*
from scipy import *

def sousblock(n,A):
    """A=np.array([[2,1,1],[1,2,1],[1,1,2]])"""
    a=np.zeros((n,n))
    for k in range(n-1):
        for i in range(2):
            for j in range(2):
                a[i+k,j+k]=a[i+k,j+k]+A[i,j]
    return(a)


# C=soublock(6,np.array([[4000/3,-5000/6],[-5000/6,4000/3]]))

#Paramètres
r=0.5e-3
H=5e-3
hi=0.5e-3

h=100
P=2*np.pi*r
Lambda=40
S=np.pi*r**2

T0=100
Tinf=20

m=h*P/(Lambda*S)

K=1/(3*hi)*np.array([[7,-8,1],[-8,16,-8],[1,-8,7]])
M=m*hi/30*np.array([[4,2,-1],[2,16,2],[-1,2,4]])

A=sousblock(11,K+M)#Matrice global
A[0,0]=1
A[0,1]=0
A[0,2]=0

B=np.array([[100],[4/6*hi*m*Tinf],[2/6*hi*m*Tinf],[4/6*hi*m*Tinf],[2/6*hi*m*Tinf],[4/6*hi*m*Tinf],[2/6*hi*m*Tinf],[4/6*hi*m*Tinf],[2/6*hi*m*Tinf],[4/6*hi*m*Tinf],[hi*m*Tinf]])

#Résolution
T=np.linalg.inv(A).dot(B)

def Tex(x):
    m=np.sqrt(h*P/(Lambda*S))
    y=Tinf+(T0-Tinf)*np.cosh(m*(H-x))/np.cosh(m*H)
    return(y)

#Tracé de la température
x=np.array([[0*hi],[1*hi],[2*hi],[3*hi],[4*hi],[5*hi],[6*hi],[7*hi],[8*hi],[9*hi],[10*hi]])
plt.plot(x,T,label='T')
plt.plot(x,Tex(x),label='Tex')
plt.title('Problème de conduction convection avec condition adiabatique')
plt.xlabel('x')
plt.ylabel('T(x)')
plt.grid(True,which="both", linestyle='--')
plt.legend()
plt.show() # affiche la figure a l'ecran

#Tracé de l'erreur l'erreur relative epsilon
Eps=np.zeros(5)
for i in range(5):
    Eps=abs(T-Tex(x))/abs(Tex(x))
plt.plot(x,Eps,label="err")
plt.loglog(x,x*x,label="pente 2")
plt.title("Courbe de l'erreur relative")
plt.xlabel('x')
plt.ylabel('Eps(x)')
plt.xscale('log')
plt.legend()
plt.grid(True,which="both", linestyle='--')
plt.show() # affiche la figure a l'ecran