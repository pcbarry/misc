#!/usr/bin/env python
import sys,os
import numpy as np
import mpmath
from scipy.special import gamma

from mpmath import fp
from qcdlib import alphaS,aux
from tools.config import conf

class DMELLIN:

    def __init__(self,muF2vals,nptsN=8,nptsM=8,extN=True,extM=True):

        xN,wN=np.polynomial.legendre.leggauss(nptsN)
        xM,wM=np.polynomial.legendre.leggauss(nptsM)
        znodesN=[0,0.1,0.3,0.6,1.0,1.6,2.4,3.5,5,7,10,14,19,25,32,40,50,63]
        znodesM=[0,0.1,0.3,0.6,1.0,1.6,2.4,3.5,5,7,10,14,19,25,32,40,50,63]
        if extN: znodesN.extend([70,80,90,100])
        if extM: znodesM.extend([70,80,90,100])

        ZN,WN,JACN=[],[],[]
        ZM,WM,JACM=[],[],[]
        for i in range(len(znodesN)-1):
            aN,bN=znodesN[i],znodesN[i+1]
            ZN.extend(0.5*(bN-aN)*xN+0.5*(aN+bN))
            WN.extend(wN)
            JACN.extend([0.5*(bN-aN) for j in range(xN.size)])
        for i in range(len(znodesM)-1):
            aM,bM=znodesM[i],znodesM[i+1]
            ZM.extend(0.5*(bM-aM)*xM+0.5*(aM+bM))
            WM.extend(wM)
            JACM.extend([0.5*(bM-aM) for j in range(xM.size)])
        self.ZN=np.array(ZN)
        self.ZM=np.array(ZM)
        # globalize
        self.WN=np.array(WN)
        self.WM=np.array(WM)
        ZN=self.ZN
        ZM=self.ZM
        self.JACN=np.array(JACN)
        self.JACM=np.array(JACM)
        # gen mellin contour                                                                                       
        cN=1.9
        cM = cN
        phi=3.0/4.0*np.pi
        phi1=np.ones(len(ZN))*phi
        self.N=cN+ZN*np.exp(complex(0,phi))
        self.phase1=np.array([np.exp(complex(0,phi1[i])) for i in range(len(phi1))])

        #--do Landau pole and adjust for self.M
        alphaS=conf['alphaS']
        euler=fp.euler

        M={}
        phase2={}
        for muF2 in np.unique(muF2vals):
            Nf=alphaS.get_Nf(muF2)
            aS=alphaS.get_alphaS(muF2)
            CA=3.0
            b0=1.0/12.0/np.pi*(11*CA-2*Nf)
            Nlandau=1.0/self.N*np.exp(1.0/aS/b0-2*euler)
            if 'westmark' in conf and conf['westmark']==True:
                Nlandau=1.0/self.N*np.exp(1.0/aS/b0)
            thetac=-np.arctan2(Nlandau.imag,Nlandau.real-cN)
            phi2nominal=phi1
            phi2other=(thetac+np.pi)/2
            phi2=[]
            M[muF2]=[]
            for i in range(len(thetac)):
                phi2.append(max(phi2nominal[i],phi2other[i]))
                M[muF2].append(cM+ZM[i]*np.exp(complex(0,phi2[i])))
            M[muF2]=np.array(M[muF2])
            phase2[muF2]=np.array([np.exp(complex(0,phi2[i])) for i in range(len(phi2))])
        self.phase2=phase2
        self.M=M

    def invert(self,Q2,f1,f2,F,G):
        f2s=np.conj(f2)
        WN=self.WN*self.JACN
        WM=self.WM*self.JACM
        xsec1=np.einsum('k,l,k,l,k,l,kl',WN,WM,f1,f2,self.phase1,self.phase2[Q2],F,optimize=True)
        xsec2=np.einsum('k,l,k,l,k,l,kl',WN,WM,f1,f2s,self.phase1,np.conj(self.phase2[Q2]),G,optimize=True)
        return np.real(-1.0/2.0/np.pi**2*(xsec1-xsec2))

if __name__=='__main__':

    conf['order']='NLO'
    conf['Q20']=1.3**2
    conf['aux']=aux.AUX()
    conf['alphaS']=alphaS.ALPHAS()

    q2=10.0
    mell=DMELLIN(q2)
    a=-1.8
    b=6.0
    N=mell.N
    M=mell.M
    Ms=np.conj(M[q2])

    m=M[q2]
    ms=Ms
    IN=np.ones(len(N))
    IM=np.ones(len(m))
    IMs=np.ones(len(ms))
    momN=gamma(N+a)*gamma(b+1)/gamma(N+a+b+1)
    momM=gamma(m+a)*gamma(b+1)/gamma(m+a+b+1)
    momMs=gamma(ms+a)*gamma(b+1)/gamma(ms+a+b+1)
    momN=np.einsum('i,j->ij',momN,IM)
    momM=np.einsum('i,j->ij',IN,momM)
    momMs=np.einsum('i,j->ij',IN,momMs)
    X=10**np.linspace(-5,-1,10)
    Y=10**np.linspace(-5,-1,10)
    f=lambda x: x**a*(1-x)**b
    g=lambda y: y**a*(1-y)**b
    for x in X:
        for y in Y:
          #print 'x=%10.4e  y=%10.4e   f=%10.4e   g=%10.4e  fg=%10.4e  inv=%10.4e'%(x,y,f(x),g(y),f(x)*g(y),mell.invert(x,y,momN*momM,momN*momMs))
          inv=mell.invert(x**(-N),y**(-m),momN*momM,momN*momMs)
          print 'x=%10.4e  y=%10.4e   f=%10.4e   g=%10.4e  fg=%10.4e  inv=%10.4e  rat=%10.4e'%(x,y,f(x),g(y),f(x)*g(y),inv,f(x)*g(y)/inv)

















