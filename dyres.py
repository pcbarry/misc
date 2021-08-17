#!/usr/bin/env python
import numpy as np
#import pylab as py
import math

from mpmath           import fp
from scipy.special    import gamma as Gamma
from scipy.integrate  import quad
from scipy.integrate  import fixed_quad
from tools.tools      import load,save,checkdir#,load_config

#from scipy.interpolate import CubicSpline as inter
from scipy.interpolate import griddata

from qcdlib import alphaS,aux,eweak
from qcdlib import mellin
from dyresmellin import DYMELLIN
from qcdlib.pdfpionalt import PDF
#from qcdlib.npdfalt import PDF as NPDF
from qcdlib.wpdfalt import PDF as WPDF
from tools.config import conf,load_config

class DYRES:

    def __init__(self):

        self.setup_constants()
        self.storage={}
        self.data={}
        self.storage_match={}
        self.alphaS=conf['alphaS']

        self.mell={}
        self.mell[+1]=conf['dymell1']
        self.mell[+2]=conf['dymell2']
        self.mell[+3]=conf['dymell3']

        #self.pdfA={}
        #self.pdfA[+1]=conf['resummed pdf-pion 1']
        #self.pdfA[+2]=conf['resummed pdf-pion 2']
        #self.pdfA[+3]=conf['resummed pdf-pion 3']
        self.pdfA=conf['pdf-pion']

        #self.pdfB={}
        #self.pdfB[-1]=conf['resummed W pdf 1']
        #self.pdfB[-2]=conf['resummed W pdf 2']
        #self.pdfB[-3]=conf['resummed W pdf 3']
        self.pdfB=conf['pdf-W']

    def setup_constants(self):
        self.CA=3.0
        self.CF=4.0/3.0
        self.euler=fp.euler

        self.eU2=4.0/9.0
        self.eD2=1.0/9.0

    def get_lum(self,fA,fB,ic,muR2,muF2):
        Nf=self.alphaS.get_Nf(muR2)

        out= self.eU2*(fA.res_storage[(ic,muF2)]['u']  * fB.res_storage[(-ic,muF2)]['ub']\
                      +fA.res_storage[(ic,muF2)]['ub'] * fB.res_storage[(-ic,muF2)]['u'])\
            +self.eD2*(fA.res_storage[(ic,muF2)]['d']  * fB.res_storage[(-ic,muF2)]['db']\
                      +fA.res_storage[(ic,muF2)]['db'] * fB.res_storage[(-ic,muF2)]['d'])\
            +self.eD2*(fA.res_storage[(ic,muF2)]['s']  * fB.res_storage[(-ic,muF2)]['sb']\
                      +fA.res_storage[(ic,muF2)]['sb'] * fB.res_storage[(-ic,muF2)]['s'])
        if Nf>3: out+=self.eU2*(fA.res_storage[(ic,muF2)]['c']  * fB.res_storage[(-ic,muF2)]['cb']\
                      +fA.res_storage[(ic,muF2)]['cb'] * fB.res_storage[(-ic,muF2)]['c'])
        if Nf>4: out+=self.eD2*(fA.res_storage[(ic,muF2)]['b']  * fB.res_storage[(-ic,muF2)]['bb']\
                      +fA.res_storage[(ic,muF2)]['bb'] * fB.res_storage[(-ic,muF2)]['b'])

        return out

    def _get_omega(self,N,muR2,muF2):
        CA,CF,euler=self.CA,self.CF,self.euler

        Nf=self.alphaS.get_Nf(muR2)
        aS=self.alphaS.get_alphaS(muR2)
        Nbar=N*np.exp(euler)

        b0=1.0/12.0/np.pi*(11*CA-2*Nf)
        b1=1.0/24.0/np.pi**2*(17*CA**2-5*CA*Nf-3*CF*Nf)
        lam=b0*aS*np.log(Nbar)
        Aq1=CF
        Aq2=1.0/2.0*CF*(CA*(67.0/18.0-np.pi**2/6.0)-5.0/9.0*Nf)

        Cq=aS/np.pi*CF*(-4+2*np.pi**2/3.0)
        h1=Aq1/2.0/np.pi/b0/lam*(2*lam+(1-2*lam)*np.log(1-2*lam))
        h2=(np.pi*Aq1*b1-b0*Aq2)*(2*lam+np.log(1-2*lam))/2.0/np.pi**2/b0**3\
           +Aq1*b1/4.0/np.pi/b0**3*np.log(1-2*lam)**2

        #match=1.0+aS/6.0/np.pi*(CF*(-24+4*np.pi**2)+12*Aq1*np.log(Nbar)**2)

        if muF2!=muR2:
            Cq+=aS/np.pi*CF*3.0/2.0*np.log(muR2/muF2)
            h2+=Aq1/2.0/np.pi/b0*np.log(1-2*lam)*np.log(muR2/muF2)
            #match+=aS/np.pi*(CF*3.0/2.0*np.log(muR2/muF2)-2*Aq1*np.log(Nbar)*np.log(muR2/muF2))

        if conf['Cq']=='exponential':
            logomega=Cq+2*h1*np.log(Nbar)+2*h2
            return np.exp(logomega)#-match
        elif conf['Cq']=='prefactor':
            logomega=2*h1*np.log(Nbar)+2*h2
            return (1+Cq)*np.exp(logomega)#-match

    def get_omega(self,muR2,muF2):
        if (muR2,muF2) not in self.storage: 
            self.storage[(muR2,muF2)]={}
            if conf['MP method']=='cosine':
              self.storage[(muR2,muF2)][1]=0.5*(self._get_omega(self.mell[1].Np,muR2,muF2)\
                                               +self._get_omega(self.mell[1].Nm,muR2,muF2))
              self.storage[(muR2,muF2)][2]=0.5*(self._get_omega(self.mell[2].Np,muR2,muF2)\
                                               +self._get_omega(self.mell[2].Nm,muR2,muF2))
              self.storage[(muR2,muF2)][3]=0.5*(self._get_omega(self.mell[3].Np,muR2,muF2)\
                                               +self._get_omega(self.mell[3].Nm,muR2,muF2))
            elif conf['MP method']=='expansion':
              self.storage[(muR2,muF2)][1]=self._get_omega(self.mell[1].N,muR2,muF2)
              self.storage[(muR2,muF2)][2]=self._get_omega(self.mell[2].N,muR2,muF2)
              self.storage[(muR2,muF2)][3]=self._get_omega(self.mell[3].N,muR2,muF2)
            elif conf['MP method']=='exact':
              self.storage[(muR2,muF2)][1]=self._get_omega_exact(self.mell[1].Np,self.mell[1].Nm,muR2,muF2)
              self.storage[(muR2,muF2)][2]=self._get_omega_exact(self.mell[2].Np,self.mell[2].Nm,muR2,muF2)
              self.storage[(muR2,muF2)][3]=self._get_omega_exact(self.mell[3].Np,self.mell[3].Nm,muR2,muF2)
        return self.storage[(muR2,muF2)]

    def _get_match(self,N,muR2,muF2):
        CA,CF,euler=self.CA,self.CF,self.euler

        Nf=self.alphaS.get_Nf(muR2)
        aS=self.alphaS.get_alphaS(muR2)
        Nbar=N*np.exp(euler)

        Aq1=CF
        match=1.0+aS/6.0/np.pi*(CF*(-24+4*np.pi**2)+12*Aq1*np.log(Nbar)**2)

        if muF2!=muR2:
            match+=aS/np.pi*(CF*3.0/2.0*np.log(muR2/muF2)-2*Aq1*np.log(Nbar)*np.log(muR2/muF2))

        return match

    def get_match(self,muR2,muF2):
        if (muR2,muF2) not in self.storage_match:
            self.storage_match[(muR2,muF2)]={}
            if conf['MP method']=='cosine':
              self.storage_match[(muR2,muF2)][1]=0.5*(self._get_match(self.mell[1].Np,muR2,muF2)\
                                                     +self._get_match(self.mell[1].Nm,muR2,muF2))
              self.storage_match[(muR2,muF2)][2]=0.5*(self._get_match(self.mell[2].Np,muR2,muF2)\
                                                     +self._get_match(self.mell[2].Nm,muR2,muF2))
              self.storage_match[(muR2,muF2)][3]=0.5*(self._get_match(self.mell[3].Np,muR2,muF2)\
                                                     +self._get_match(self.mell[3].Nm,muR2,muF2))
            elif conf['MP method']=='expansion':
              self.storage_match[(muR2,muF2)][1]=self._get_match(self.mell[1].N,muR2,muF2)
              self.storage_match[(muR2,muF2)][2]=self._get_match(self.mell[2].N,muR2,muF2)
              self.storage_match[(muR2,muF2)][3]=self._get_match(self.mell[3].N,muR2,muF2)
        return self.storage_match[(muR2,muF2)]

    def _get_xsec(self,tau,muR2,muF2,savedir=False,savename=None,verb=False):
        if verb: print 'computing omega'
        omega=self.get_omega(muR2,muF2)
        match=self.get_match(muR2,muF2)
        if verb: print 'evolving pdfs'
        #for _ in self.pdfA:  self.pdfA[_].evolve(muF2) 
        #for _ in self.pdfB:  self.pdfB[_].evolve(muF2)
        self.pdfA.evolve(muF2)
        self.pdfB.evolve(muF2)
        if savedir: 
            data=self.data
            data['muF2=%.2f'%muF2]={}
            data['muF2=%.2f'%muF2]['NM']={}
            data['muF2=%.2f'%muF2]['M']={}
            data['muF2=%.2f'%muF2]['match']={}


        if verb: print 'integrating over N'
        sigM=0
        matchM=0
        for ic in [1,2,3]:
            lum=self.get_lum(self.pdfA,self.pdfB,ic,muR2,muF2)

            sigNM = lum*omega[ic]
            matchNM=lum*match[ic]
            if ic==1:
                sigM-=self.mell[ic].invert(tau,sigNM)
                matchM-=self.mell[ic].invert(tau,matchNM)
                if savedir:
                    data['muF2=%.2f'%muF2]['NM'][ic]=-sigNM
                    data['muF2=%.2f'%muF2]['M'][ic]=-self.mell[ic].invert(tau,sigNM)
                    data['muF2=%.2f'%muF2]['match'][ic]=-self.mell[ic].invert(tau,matchNM)

            else:
                sigM+=self.mell[ic].invert(tau,sigNM)
                matchM+=self.mell[ic].invert(tau,matchNM)
                if savedir:
                    data['muF2=%.2f'%muF2]['NM'][ic]=sigNM
                    data['muF2=%.2f'%muF2]['M'][ic]=self.mell[ic].invert(tau,sigNM)
                    data['muF2=%.2f'%muF2]['match'][ic]=self.mell[ic].invert(tau,matchNM)

        totsig=sigM-matchM

        if verb: print 'integrating over M'
        #return sig
        if savedir:
            self.data=data
            checkdir('%s'%savedir)
            save(self.data,'%s/%s'%(savedir,savename))
        return totsig

    def interpolating_function(self,tau,muR2,muF2,M):

      aEM=conf['eweak'].get_alpha(muR2)
      S=muF2/tau
      prefactor=4*np.pi*aEM**2/9.0/S/muF2

      xsec=self._get_xsec(tau,muR2,muF2)
      raw=prefactor*xsec

      return griddata(conf['dymell1'].M,raw,M,fill_value=0,method='cubic')

    def get_xsec(self,S,Y,muR2,muF2):
      tau=muF2/S
      func=lambda M: np.exp(-1j*M*Y)/2./np.pi*self.interpolating_function(tau,muR2,muF2,M)
      return 2*np.real(fixed_quad(func,0,30,n=200)[0])


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

    conf['dymell1']=DYMELLIN(contour=1)
    conf['dymell2']=DYMELLIN(contour=2)
    conf['dymell3']=DYMELLIN(contour=3)

    conf['pdf-pion']=PDF()

    conf['pdf-W']=WPDF()

    dyres=DYRES()

    tau=0.3345
    muR2=53.05
    muF2=53.05
    S=muF2/tau
    Y=0.074669
    print dyres.get_interpolated(S,Y,muR2,muF2)
    print dyres.get_interpolated(S,0.1,muR2,muF2)
    #dyres.get_xsec(tau,muR2,muF2,verb=True)




