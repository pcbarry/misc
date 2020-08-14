#!/usr/bin/env python
import numpy as np
#import pylab as py
import math

from mpmath           import fp
from scipy.special    import gamma as Gamma
from scipy.integrate  import quad
from scipy.integrate  import fixed_quad
from tools.tools      import load,save#,checkdir,load_config

from scipy.interpolate import CubicSpline as inter

from qcdlib import alphaS,aux,eweak
from qcdlib import mellin
from dyresmellinexact import DMELLIN
from qcdlib.pdfpionaltexact import PDF
from qcdlib.wpdfaltexact import PDF as WPDF
from tools.config import conf

class DYRESEXACT:

    def __init__(self):

        self.setup_constants()
        self.storage={}
        self.alphaS=conf['alphaS']

        self.mell=conf['double mellin resum']

        self.pdfA=conf['pdf-pion']

        self.pdfB=conf['pdf-W']

    def setup_constants(self):
        self.CA=3.0
        self.CF=4.0/3.0
        self.euler=fp.euler

        self.eU2=4.0/9.0
        self.eD2=1.0/9.0

    def get_lum(self,fA,fB,muR2,muF2,Mconj=False):
        Nf=self.alphaS.get_Nf(muR2)

       # if Mconj==True:
       #     out= self.eU2*(fA.storage[muF2]['u']  * np.conj(fB.res_storage[muF2]['ub'])\
       #                   +fA.storage[muF2]['ub'] * np.conj(fB.res_storage[muF2]['u']))\
       #         +self.eD2*(fA.storage[muF2]['d']  * np.conj(fB.res_storage[muF2]['db'])\
       #                   +fA.storage[muF2]['db'] * np.conj(fB.res_storage[muF2]['d']))\
       #         +self.eD2*(fA.storage[muF2]['s']  * np.conj(fB.res_storage[muF2]['sb'])\
       #                   +fA.storage[muF2]['sb'] * np.conj(fB.res_storage[muF2]['s']))
       #     if Nf>3: out+=self.eU2*(fA.storage[muF2]['c']  * np.conj(fB.res_storage[muF2]['cb'])\
       #                   +fA.storage[muF2]['cb'] * np.conj(fB.res_storage[muF2]['c']))
       #     if Nf>4: out+=self.eD2*(fA.storage[muF2]['b']  * np.conj(fB.res_storage[muF2]['bb'])\
       #                   +fA.storage[muF2]['bb'] * np.conj(fB.res_storage[muF2]['b']))

       # elif Mconj==False:
       #     out= self.eU2*(fA.storage[muF2]['u']  * fB.res_storage[muF2]['ub']\
       #                   +fA.storage[muF2]['ub'] * fB.res_storage[muF2]['u'])\
       #         +self.eD2*(fA.storage[muF2]['d']  * fB.res_storage[muF2]['db']\
       #                   +fA.storage[muF2]['db'] * fB.res_storage[muF2]['d'])\
       #         +self.eD2*(fA.storage[muF2]['s']  * fB.res_storage[muF2]['sb']\
       #                   +fA.storage[muF2]['sb'] * fB.res_storage[muF2]['s'])
       #     if Nf>3: out+=self.eU2*(fA.storage[muF2]['c']  * fB.res_storage[muF2]['cb']\
       #                   +fA.storage[muF2]['cb'] * fB.res_storage[muF2]['c'])
       #     if Nf>4: out+=self.eD2*(fA.storage[muF2]['b']  * fB.res_storage[muF2]['bb']\
       #                   +fA.storage[muF2]['bb'] * fB.res_storage[muF2]['b'])

        if Mconj==True:
            out= self.eU2*(np.einsum('i,j->ij',fA.storage[muF2]['u'],np.conj(fB.res_storage[muF2]['ub']))\
                          +np.einsum('i,j->ij',fA.storage[muF2]['ub'],np.conj(fB.res_storage[muF2]['u'])))\
                +self.eD2*(np.einsum('i,j->ij',fA.storage[muF2]['d'],np.conj(fB.res_storage[muF2]['db']))\
                          +np.einsum('i,j->ij',fA.storage[muF2]['db'],np.conj(fB.res_storage[muF2]['d'])))\
                +self.eD2*(np.einsum('i,j->ij',fA.storage[muF2]['s'],np.conj(fB.res_storage[muF2]['sb']))\
                          +np.einsum('i,j->ij',fA.storage[muF2]['sb'],np.conj(fB.res_storage[muF2]['s'])))
            if Nf>3: out+=self.eU2*(np.einsum('i,j->ij',fA.storage[muF2]['c'],np.conj(fB.res_storage[muF2]['cb']))\
                          +np.einsum('i,j->ij',fA.storage[muF2]['cb'],np.conj(fB.res_storage[muF2]['c'])))
            if Nf>4: out+=self.eD2*(np.einsum('i,j->ij',fA.storage[muF2]['b'],np.conj(fB.res_storage[muF2]['bb']))\
                          +np.einsum('i,j->ij',fA.storage[muF2]['bb'],np.conj(fB.res_storage[muF2]['b'])))

        elif Mconj==False:
            out= self.eU2*(np.einsum('i,j->ij',fA.storage[muF2]['u'],fB.res_storage[muF2]['ub'])\
                          +np.einsum('i,j->ij',fA.storage[muF2]['ub'],fB.res_storage[muF2]['u']))\
                +self.eD2*(np.einsum('i,j->ij',fA.storage[muF2]['d'],fB.res_storage[muF2]['db'])\
                          +np.einsum('i,j->ij',fA.storage[muF2]['db'],fB.res_storage[muF2]['d']))\
                +self.eD2*(np.einsum('i,j->ij',fA.storage[muF2]['s'],fB.res_storage[muF2]['sb'])\
                          +np.einsum('i,j->ij',fA.storage[muF2]['sb'],fB.res_storage[muF2]['s']))
            if Nf>3: out+=self.eU2*(np.einsum('i,j->ij',fA.storage[muF2]['c'],fB.res_storage[muF2]['cb'])\
                          +np.einsum('i,j->ij',fA.storage[muF2]['cb'],fB.res_storage[muF2]['c']))
            if Nf>4: out+=self.eD2*(np.einsum('i,j->ij',fA.storage[muF2]['b'],fB.res_storage[muF2]['bb'])\
                          +np.einsum('i,j->ij',fA.storage[muF2]['bb'],fB.res_storage[muF2]['b']))
        return out

    def _get_omega(self,N1,N2,muR2,muF2):
        CA,CF,euler=self.CA,self.CF,self.euler

        Nf=self.alphaS.get_Nf(muR2)
        aS=self.alphaS.get_alphaS(muR2)
        N1bar=N1*np.exp(euler)
        N2bar=N2*np.exp(euler)

        b0=1.0/12.0/np.pi*(11*CA-2*Nf)
        b1=1.0/24.0/np.pi**2*(17*CA**2-5*CA*Nf-3*CF*Nf)
        I1=np.ones(len(N1bar))
        I2=np.ones(len(N2bar))
        #lam1=b0*aS*np.log(N1bar)
        #lam2=b0*aS*np.log(N2bar)
        lam1=b0*aS*np.einsum('i,j->ij',np.log(N1bar),I2)
        lam2=b0*aS*np.einsum('i,j->ij',I1,np.log(N2bar))
        Aq1=CF
        Aq2=1.0/2.0*CF*(CA*(67.0/18.0-np.pi**2/6.0)-5.0/9.0*Nf)

        Cq=aS/np.pi*CF*(-4+2*np.pi**2/3.0)
        h1=Aq1/np.pi/b0/(lam1+lam2)*(lam1+lam2+(1-lam1-lam2)*np.log(1-lam1-lam2))
        h2=Aq1*b1/2.0/np.pi/b0**3*(2*(lam1+lam2+lam1*np.log(1-lam1)+lam2*np.log(1-lam2))\
                                   +(2-lam1)*np.log(1-lam1)**2+(2-lam2)*np.log(1-lam2)**2\
                                   +(1-lam1-lam2)*(2-np.log(1-lam1-lam2))*np.log(1-lam1-lam2))\
                -Aq2/b0**2/np.pi**2*(lam1+lam2+np.log(1-lam1-lam2))

        match=1.0+Aq1*aS/6.0/np.pi*(-24+4*np.pi**2+3*np.log(N1bar)**2+3*np.log(N2bar)**2+6*np.log(N1bar)*np.log(N2bar))

        if muF2!=muR2:
            Cq+=aS/np.pi*CF*3.0/2.0*np.log(muR2/muF2)
            h2+=Aq1/b0/np.pi*(np.log(1-lam1-lam2)*np.log(muR2/muF2))
            match+=Aq1*aS/6.0/np.pi*(9*np.log(muR2/muF2)-6*(np.log(N1bar)+np.log(N2bar))*np.log(muR2/muF2))

        logomega=Cq+h1*(np.log(N1bar)+np.log(N2bar))+h2
        return np.exp(logomega)-match

    def get_omega(self,muR2,muF2):
        if (muR2,muF2) not in self.storage: 
            self.storage[(muR2,muF2)]={}
            self.storage[(muR2,muF2)]['NM']=self._get_omega(self.mell.N,self.mell.M[muF2],muR2,muF2)
            self.storage[(muR2,muF2)]['NMs']=self._get_omega(self.mell.N,np.conj(self.mell.M[muF2]),muR2,muF2)
        return self.storage[(muR2,muF2)]

    def get_xsec(self,tau,Y,muR2,muF2,verb=False):
        if verb: print 'computing omega'
        omega=self.get_omega(muR2,muF2)
        if verb: print 'evolving pdfs'
        self.pdfA.evolve(muF2)
        self.pdfB.evolve(muF2)

        lumNM=self.get_lum(self.pdfA,self.pdfB,muR2,muF2,Mconj=False)
        lumNMs=self.get_lum(self.pdfA,self.pdfB,muR2,muF2,Mconj=True)

        sigNM = lumNM*omega['NM']
        sigNMs= lumNMs*omega['NMs']

        if verb: print 'inverting'
        x10term=(tau**0.5*np.exp(Y))**(-self.mell.N)
        x20term=(tau**0.5*np.exp(-Y))**(-self.mell.M[muF2])
        sig=self.mell.invert(muR2,x10term,x20term,sigNM,sigNMs)

        aEM=conf['eweak'].get_alpha(muR2)
        S=muF2/tau
        prefactor=4*np.pi*aEM**2/9.0/S/muF2
        xsec=prefactor*sig
        return xsec

if __name__=="__main__":

    conf['order']='NLO'
    conf['Q20'] = 1.0 
    conf['aux']=aux.AUX()
    conf['dglap mode']='truncated'
    conf['alphaS']=alphaS.ALPHAS()
    conf['Z']=74.
    conf['A']=184.
    conf['eweak']=eweak.EWEAK()
    conf['mellin-pion']=mellin.MELLIN(npts=8,extended=True)
    conf['mellin']=mellin.MELLIN()
    conf['MP method']='cosine'

    tau=0.3345**2
    muR2=53.05
    muF2=53.05
    S=muF2/tau
    Y=0.074669

    conf['double mellin resum']=DMELLIN(np.array(muF2))

    conf['pdf-pion']=PDF()

    conf['pdf-W']=WPDF()

    dyres=DYRESEXACT()

    print dyres.get_xsec(tau,Y,muR2,muF2,verb=True)
    print dyres.get_xsec(tau,0.1,muR2,muF2,verb=True)
    #dyres.get_xsec(tau,muR2,muF2,verb=True)




