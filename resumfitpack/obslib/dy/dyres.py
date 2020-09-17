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
from dyresmellin import DYMELLIN
from qcdlib.pdfpionalt import PDF
#from qcdlib.npdfalt import PDF as NPDF
from qcdlib.wpdfalt import PDF as WPDF
from tools.config import conf

class DYRES:

    #@profile
    def __init__(self):

        self.setup_constants()
        self.storage={}
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

    #@profile
    def setup_constants(self):
        self.CA=3.0
        self.CF=4.0/3.0
        self.euler=fp.euler

        self.eU2=4.0/9.0
        self.eD2=1.0/9.0

    #@profile
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

    def _get_westmark_omega(self,N,muR2,muF2):
        CA,CF,euler=self.CA,self.CF,self.euler

        Nf=self.alphaS.get_Nf(muR2)
        aS=self.alphaS.get_alphaS(muR2)

        b0=1.0/12.0/np.pi*(11*CA-2*Nf)
        b1=1.0/24.0/np.pi**2*(17*CA**2-5*CA*Nf-3*CF*Nf)/b0
        lam=b0*aS*np.log(N)
        Aq1=CF/np.pi
        Aq2=1.0/2.0*CF*(CA*(67.0/18.0-np.pi**2/6.0)-5.0/9.0*Nf)/np.pi**2

        #Cq=aS/np.pi*CF*(-4+2*np.pi**2/3.0)
        g01=CF/2.0/np.pi*(4*np.pi**2/3-8+4*euler**2)
        g1=Aq1/b0/lam*(2*lam+(1-2*lam)*np.log(1-2*lam))
        g2=Aq1/b0*(np.log(1-2*lam)*(-2*euler))\
                +Aq1*b1/b0**2*(0.5*np.log(1-2*lam)**2+np.log(1-2*lam)+2*lam)\
                -Aq2/b0**2*(2*lam+np.log(1-2*lam))

        match=1.0+aS*(g01+Aq1*(2*np.log(N)**2+2*np.log(N)*2*euler))

        if muF2!=muR2:
            g01+=-2*CF/np.pi*euler*np.log(muR2/muF2)
            h2+=Aq1/b0*np.log(1-2*lam)*np.log(muR2/muF2)
            match+=-aS*Aq1*2*np.log(N)*np.log(muR2/muF2)

        g0=1+aS*g01
        logomega=np.log(N)*g1+g2
        return g0*np.exp(logomega)-match

    #@profile
    def _get_omega(self,N,muR2,muF2):
        if 'westmark' in conf and conf['westmark']==True: return self._get_westmark_omega(N,muR2,muF2)
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

        match=1.0+aS/6.0/np.pi*(CF*(-24+4*np.pi**2)+12*Aq1*np.log(Nbar)**2)

        if muF2!=muR2:
            Cq+=aS/np.pi*CF*3.0/2.0*np.log(muR2/muF2)
            h2+=Aq1/2.0/np.pi/b0*np.log(1-2*lam)*np.log(muR2/muF2)
            match+=aS/np.pi*(CF*3.0/2.0*np.log(muR2/muF2)-2*Aq1*np.log(Nbar)*np.log(muR2/muF2))

        if conf['Cq']=='exponential':
            logomega=Cq+2*h1*np.log(Nbar)+2*h2
            return np.exp(logomega)-match
        elif conf['Cq']=='prefactor':
            logomega=2*h1*np.log(Nbar)+2*h2
            return (1+Cq)*np.exp(logomega)-match

    def _get_omega_exact(self,N1,N2,muR2,muF2):
        CA,CF,euler=self.CA,self.CF,self.euler

        Nf=self.alphaS.get_Nf(muR2)
        aS=self.alphaS.get_alphaS(muR2)
        N1bar=N1*np.exp(euler)
        N2bar=N2*np.exp(euler)

        b0=1.0/12.0/np.pi*(11*CA-2*Nf)
        b1=1.0/24.0/np.pi**2*(17*CA**2-5*CA*Nf-3*CF*Nf)
        lam1=b0*aS*np.log(N1bar)
        lam2=b0*aS*np.log(N2bar)
        Aq1=CF
        Aq2=1.0/2.0*CF*(CA*(67.0/18.0-np.pi**2/6.0)-5.0/9.0*Nf)

        Cq=aS/np.pi*CF*(-4+2*np.pi**2/3.0)
        #h1=Aq1/2.0/np.pi/b0/lam*(lam+(1-lam)*np.log(1-lam))
        h1=Aq1/np.pi/b0/(lam1+lam2)*(lam1+lam2+(1-lam1-lam2)*np.log(1-lam1-lam2))
        #h2=1.0/2.0/np.pi**2/b0**3*(np.pi*Aq1*b1*(-lam+np.log(1-lam))-b0*Aq2*(lam+np.log(1-lam)))\
        #h2=1.0/2.0/np.pi**2/b0**3*(np.pi*Aq1*b1*(lam+np.log(1-lam))-b0*Aq2*(lam+np.log(1-lam)))\
        #   +Aq1*b1/4.0/np.pi/b0**3*np.log(1-lam)**2
        h2=Aq1*b1/2.0/np.pi/b0**3*(2*(lam1+lam2+lam1*np.log(1-lam1)+lam2*np.log(1-lam2))\
                                   +(2-lam1)*np.log(1-lam1)**2+(2-lam2)*np.log(1-lam2)**2\
                                   +(1-lam1-lam2)*(2-np.log(1-lam1-lam2))*np.log(1-lam1-lam2))\
                -Aq2/b0**2/np.pi**2*(lam1+lam2+np.log(1-lam1-lam2))

        #match=1.0+aS/6.0/np.pi*(CF*(-24+4*np.pi**2)+3*Aq1*(np.log(N1bar)+np.log(N2bar))**2\
        #                       -12*Aq1*b1/b0**2*(np.log(N1bar)+np.log(N2bar)))
        match=1.0+Aq1*aS/6.0/np.pi*(-24+4*np.pi**2+3*np.log(N1bar)**2+3*np.log(N2bar)**2+6*np.log(N1bar)*np.log(N2bar))

        if muF2!=muR2:
            Cq+=aS/np.pi*CF*3.0/2.0*np.log(muR2/muF2)
            #h2+=Aq1/2.0/np.pi/b0*np.log(1-lam)*np.log(muR2/muF2)
            h2+=Aq1/b0/np.pi*(np.log(1-lam1-lam2)*np.log(muR2/muF2))
            #match+=aS/np.pi*(CF*3.0/2.0*np.log(muR2/muF2)-Aq1*(np.log(N1bar)+np.log(N2bar))*np.log(muR2/muF2))
            match+=Aq1*aS/6.0/np.pi*(9*np.log(muR2/muF2)-6*(np.log(N1bar)+np.log(N2bar))*np.log(muR2/muF2))

        logomega=Cq+h1*(np.log(N1bar)+np.log(N2bar))+h2
        return np.exp(logomega)-match

    #@profile
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

    #@profile
    def get_xsec(self,tau,muR2,muF2,verb=False):
        if verb: print 'computing omega'
        omega=self.get_omega(muR2,muF2)
        if verb: print 'evolving pdfs'
        #for _ in self.pdfA:  self.pdfA[_].evolve(muF2) 
        #for _ in self.pdfB:  self.pdfB[_].evolve(muF2)
        self.pdfA.evolve(muF2)
        self.pdfB.evolve(muF2)

        if verb: print 'integrating over N'
        sigM=0
        for ic in [1,2,3]:
            lum=self.get_lum(self.pdfA,self.pdfB,ic,muR2,muF2)

            sigNM = lum*omega[ic]
            if ic==1: sigM-=self.mell[ic].invert(tau,sigNM)

            else: sigM+=self.mell[ic].invert(tau,sigNM)

        if verb: print 'integrating over M'
        sig=self.mell[1].fft(sigM)/(2.0*np.pi)**2
        return sig

    #@profile
    def gather_points(self,tau,muR2,muF2):

      aEM=conf['eweak'].get_alpha(muR2)
      S=muF2/tau
      prefactor=4*np.pi*aEM**2/9.0/S/muF2

      #if (muR2,muF2) not in self.storage:
      xsec=self.get_xsec(tau,muR2,muF2)
      #self.storage[(muR2,muF2)]['fft']={}
      #self.storage[(muR2,muF2)]['fft']['raw']=prefactor*xsec
      raw=prefactor*xsec

      Y=np.linspace(0,1.5,conf['dymell1'].nY)*(2*np.pi)**2
      #--Gather small Y even points      
      #self.storage[(muR2,muF2)]['fft']['even points']=[]
      even_points=[]
      #self.storage[(muR2,muF2)]['fft']['Y even']=[]
      Y_even=[]
      for i in range(len(Y)):
        if Y[i]<2.5:
          if i%2==0:
            #self.storage[(muR2,muF2)]['fft']['even points'].append(self.storage[(muR2,muF2)]['fft']['raw'][i])
            even_points.append(raw[i])
            #self.storage[(muR2,muF2)]['fft']['Y even'].append(Y[i])
            Y_even.append(Y[i])


      #--Make into np.array
      #self.storage[(muR2,muF2)]['fft']['even points']=np.array(self.storage[(muR2,muF2)]['fft']['even points'])
      even_points=np.array(even_points)
      #self.storage[(muR2,muF2)]['fft']['Y even']=np.array(self.storage[(muR2,muF2)]['fft']['Y even'])
      Y_even=np.array(Y_even)
      return Y_even,even_points

    #@profile
    def get_interpolated(self,S,Y,muR2,muF2):
      tau=muF2/S
      x,y=self.gather_points(tau,muR2,muF2)
      #interpolated=inter(self.storage[(muR2,muF2)]['fft']['Y even'],self.storage[(muR2,muF2)]['fft']['even points'])
      if len(x)<2: return 0
      for i in range(len(y)):
          if math.isnan(y[i]): y[i]=1e100
      interpolated=inter(x,y)
      return interpolated(Y).real

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




