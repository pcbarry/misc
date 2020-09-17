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
        Ymax=1.5    #--max rapidity covered by experiments
        #Ymax=1.1057    #--max rapidity covered by experiments
        #nY = 1024   #--number of fft points to be generated in rap space
        nY = 512   #--number of fft points to be generated in rap space
        nM = 512    #--number of M points to be use
        #nM = 1024    #--number of M points to be use

        #--construct rap array with linear spacing 
        Y    = np.linspace(0,Ymax,nY)
        dY   = Y[1]-Y[0]
        Mmax = 1/dY/(2*np.pi)/2 #--need to clarify where this comes from 
        self.M=np.linspace(-Mmax,Mmax,nM)
        
        #--broadcast
        self.nY=nY

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

        """
        In order to shorten the calculation of resummation,
        we exploit symmetries of the contours.
        When initializing the PDFs, one must pass through
        the N values of the contours into the kernels and dglap.

        The pdfA is evaluated at Nc+iM/2.
        When M is positive, the contour is M-indepedent,
        and we should only pass once the N1 contours for positive M.
        We will call this self.Npmin.

        The pdfB is evaluated at Nc-iM/2.
        When M is negative, the contour is M-independent,
        and we should only pass once the N1 contours for negative M.
        We will call this self.Nmmin.

        The PDFs are written to reflect this exploitation.
        """

        M=np.copy(self.M)
        negprimeM=M[:len(M)/2+1] #--all negative Ms and one positive M
        posprimeM=M[len(M)/2-1:] #--all positive Ms and one negative M

        IM=np.ones(len(M))
        InegprimeM=np.ones(len(negprimeM))
        IposprimeM=np.ones(len(posprimeM))
        IN=np.ones(len(N))
        M=np.einsum('i,j->ij',IN,M)
        negprimeM=np.einsum('i,j->ij',IN,negprimeM)
        posprimeM=np.einsum('i,j->ij',IN,posprimeM)
        negprimeN=np.einsum('i,j->ij',N,InegprimeM)
        posprimeN=np.einsum('i,j->ij',N,IposprimeM)
        N=np.einsum('i,j->ij',N,IM)

        #--bradcast shifted moments
        Nc         =  N  - 1j*np.abs(M)/2
        #self.N     = (N  - 1j*np.abs(M)/2).flatten()
        negprimeNc =  negprimeN  - 1j*np.abs(negprimeM)/2
        posprimeNc =  posprimeN  - 1j*np.abs(posprimeM)/2
        self.Np    = (Nc + 1j*M/2).flatten()
        self.Npmin = (negprimeNc + 1j*negprimeM/2).flatten()
        self.Nm    = (Nc - 1j*M/2).flatten()
        self.Nmmin = (posprimeNc - 1j*posprimeM/2).flatten()
        self.N     = Nc.flatten()

        #--build mellin inversion 
        W   = np.einsum('i,j->ij',W,IM)
        JAC = np.einsum('i,j->ij',JAC,IM)
        shape=N.shape
        #shape=self.N.shape
        #print shape
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

        """
        In order to shorten the calculation of resummation,
        we exploit symmetries of the contours.
        When initializing the PDFs, one must pass through
        the N values of the contours into the kernels and dglap.

        The pdfA is evaluated at Nc+iM/2.
        When M is negative, the contour is M-indepedent,
        and we should only pass once the N3 contours for negative M.
        We will call this self.Npmin.

        The pdfB is evaluated at Nc-iM/2.
        When M is positive, the contour is M-independent,
        and we should only pass once the N3 contours for positive M.
        We will call this self.Nmmin.

        The PDFs are written to reflect this exploitation.
        """

        M=np.copy(self.M)
        negprimeM=M[:len(M)/2+1] #--all negative Ms and one positive M
        posprimeM=M[len(M)/2-1:] #--all positive Ms and one negative M

        IM=np.ones(len(M))
        InegprimeM=np.ones(len(negprimeM))
        IposprimeM=np.ones(len(posprimeM))
        IN=np.ones(len(N))
        M=np.einsum('i,j->ij',IN,M)
        negprimeM=np.einsum('i,j->ij',IN,negprimeM)
        posprimeM=np.einsum('i,j->ij',IN,posprimeM)
        negprimeN=np.einsum('i,j->ij',N,InegprimeM)
        posprimeN=np.einsum('i,j->ij',N,IposprimeM)
        N=np.einsum('i,j->ij',N,IM)

        #--bradcast shifted moments
        #self.N=(N + sign*1j*M/2).flatten()
        Nc         =  N  + 1j*np.abs(M)/2
        self.N     = (N  + 1j*np.abs(M)/2).flatten()
        negprimeNc =  negprimeN  + 1j*np.abs(negprimeM)/2
        posprimeNc =  posprimeN  + 1j*np.abs(posprimeM)/2
        self.Np    = (Nc + 1j*M/2).flatten()
        self.Npmin = (posprimeNc + 1j*posprimeM/2).flatten()
        self.Nm    = (Nc - 1j*M/2).flatten()
        self.Nmmin = (negprimeNc - 1j*negprimeM/2).flatten()

        #--build mellin inversion 
        W   = np.einsum('i,j->ij',W,IM)
        JAC = np.einsum('i,j->ij',JAC,IM)
        shape=N.shape
        #--here F is expected to be a tensor of rank 2 (NxM)  flattend
        self.invert=lambda x,F:np.einsum('ij,ij,ij,ij->j',x**(-Nc),\
                                            F.reshape(shape),W,JAC)*phase/(2*np.pi*1j)

    #--fast fourier transform

    def _fft(self,x,n,s):
        """
        !!! NEED to improve explanation
        x : is the integrand of M
        n : ....
        s : ....

        ---old---
        n is the length of the desired X, must be an integer power of 2
        if x is of length L, s=L/n
        #if len(x)==0: pass
        """
        X=np.zeros(n,dtype=complex)
        xnew=[]
        for i in range(len(x)): 
            if i%s==0: xnew.append(x[i])
    
        x1=xnew
        x2=[]
        for j in range(len(xnew)):
            if j==0: pass
            else: x2.append(xnew[j])
        
        if n==1: X[0]=x[0]
        else:
            Xfirsthalf=self._fft(x1,n/2,2)
            Xsecondhalf=self._fft(x2,n/2,2)
            for k in range(n):
                if k<=n/2-1: 
                    X[k]=Xfirsthalf[k]
                elif k>=n/2:
                    i=k-n/2
                    X[k]=Xsecondhalf[i]
                    
            for k in range(n/2):
                t=X[k]
                X[k]=t+np.exp(-2*np.pi*1j*k/n)*X[k+n/2]
                X[k+n/2]=t-np.exp(-2*np.pi*1j*k/n)*X[k+n/2]
        return X

    def fft(self,sigM):
        return self._fft(sigM,self.nY,1)


if __name__=='__main__':


    dymellin=DYMELLIN(contour=1)

















