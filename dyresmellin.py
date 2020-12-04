#!/usr/bin/env python
import os,sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import gamma
from scipy.integrate import quad

class DYMELLIN:

    #@profile
    def __init__(self,contour,npts=4):#8):
 
        self.set_M_values() 
        if contour==1: self.set_countour1(npts)
        if contour==2: self.set_countour2(npts)
        if contour==3: self.set_countour3(npts)

    #@profile
    def set_M_values(self):
        M=np.linspace(-0.25,0.9,10)
        M=np.append(M,10**np.linspace(0,np.log10(30),40))
        self.M=M

    #@profile
    def set_countour1(self,npts):
        x,w=np.polynomial.legendre.leggauss(npts)
        znodes=[0,0.1,0.3,0.6,1.0,1.6,2.4,3.5,5,7,10,14,19,25,32,40,50,63]
        #znodes.extend([70,80,90,100])
        zlen=len(znodes)

        Z,W,JAC=[],[],[]
        for i in range(zlen-1):
            a,b=znodes[i],znodes[i+1]
            Z.extend(0.5*(b-a)*x+0.5*(a+b))
            W.extend(w)
            JAC.extend([0.5*(b-a) for j in range(x.size)])
        Z   = np.array(Z)
        W   = np.array(W)
        JAC = np.array(JAC)

        c=1.9
        phi=5*np.pi/4
        phase=np.exp(complex(0,phi))
        N=c+Z*np.exp(complex(0,phi)) 

        M=np.copy(self.M)

        IM=np.ones(len(M))
        IN=np.ones(len(N))
        M=np.einsum('i,j->ij',IN,M)
        N=np.einsum('i,j->ij',N,IM)

        #--bradcast shifted moments
        Nc         =  N  - 1j*np.abs(M)/2
        self.Np    = (Nc + 1j*M/2).flatten()
        self.Nm    = (Nc - 1j*M/2).flatten()
        self.N     = Nc.flatten()

        #--build mellin inversion 
        W   = np.einsum('i,j->ij',W,IM)
        JAC = np.einsum('i,j->ij',JAC,IM)
        shape=N.shape
        #--here F is expected to be a tensor of rank 2 (NxM)  flattend
        self.invert=lambda x,F:np.einsum('ij,ij,ij,ij->j',x**(-Nc)\
                                ,F.reshape(shape),W,JAC)*phase/(2*np.pi*1j)

    def set_countour2(self,npts):
        x,w=np.polynomial.legendre.leggauss(npts)
        znodes=np.linspace(0,1,20)
        zlen=len(znodes)

        Z,W,JAC=[],[],[]
        for i in range(zlen-1):
            a,b=znodes[i],znodes[i+1]
            Z.extend(0.5*(b-a)*x+0.5*(a+b))
            W.extend(w)
            JAC.extend([0.5*(b-a) for j in range(x.size)])
        Z   = np.array(Z)
        W   = np.array(W)
        JAC = np.array(JAC)


        M=np.copy(self.M)

        IM=np.ones(len(M))
        IZ=np.ones(len(Z))
        M=np.einsum('i,j->ij',IZ,M)
        Z=np.einsum('i,j->ij',Z,IM)

        c = 1.9
        #N = c - 1j*M/2 + Z*1j*M
        N = c - 1j*np.abs(M)/2 + Z*1j*np.abs(M)
        Nc= N

        #--bradcast shifted moments
        #self.N=(N + sign*1j*M/2).flatten()
        self.N =(N).flatten()
        self.Np=(N + 1j*M/2).flatten()
        self.Nm=(N - 1j*M/2).flatten()

        #--build mellin inversion 
        W   = np.einsum('i,j->ij',W,IM)
        JAC = np.einsum('i,j->ij',JAC,IM)
        shape=N.shape
        #--here F is expected to be a tensor of rank 2 (NxM)  flattend
        #self.invert = lambda x,F:np.einsum('ij,ij,ij,ij->j',x**(-Nc)\
        self.invert = lambda x,F:np.einsum('ij,ij,ij,ij,ij->j',x**(-Nc)\
                                   ,F.reshape(shape),W,JAC,np.abs(M))/(2*np.pi)

    def set_countour3(self,npts):
        x,w=np.polynomial.legendre.leggauss(npts)
        znodes=[0,0.1,0.3,0.6,1.0,1.6,2.4,3.5,5,7,10,14,19,25,32,40,50,63]
        #znodes.extend([70,80,90,100])
        zlen=len(znodes)

        Z,W,JAC=[],[],[]
        for i in range(zlen-1):
            a,b=znodes[i],znodes[i+1]
            Z.extend(0.5*(b-a)*x+0.5*(a+b))
            W.extend(w)
            JAC.extend([0.5*(b-a) for j in range(x.size)])
        Z   = np.array(Z)
        W   = np.array(W)
        JAC = np.array(JAC)

        c=1.9
        phi=3*np.pi/4
        phase=np.exp(complex(0,phi))
        N=c+Z*np.exp(complex(0,phi)) 

        M=np.copy(self.M)

        IM=np.ones(len(M))
        IN=np.ones(len(N))
        M=np.einsum('i,j->ij',IN,M)
        N=np.einsum('i,j->ij',N,IM)

        #--bradcast shifted moments
        Nc         =  N  + 1j*np.abs(M)/2
        self.N     = (N  + 1j*np.abs(M)/2).flatten()
        self.Np    = (Nc + 1j*M/2).flatten()
        self.Nm    = (Nc - 1j*M/2).flatten()

        #--build mellin inversion 
        W   = np.einsum('i,j->ij',W,IM)
        JAC = np.einsum('i,j->ij',JAC,IM)
        shape=N.shape
        #--here F is expected to be a tensor of rank 2 (NxM)  flattend
        self.invert=lambda x,F:np.einsum('ij,ij,ij,ij->j',x**(-Nc),\
                                            F.reshape(shape),W,JAC)*phase/(2*np.pi*1j)


if __name__=='__main__':


    dymellin=DYMELLIN(contour=1)

















