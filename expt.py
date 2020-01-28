#!/usr/bin/env python
import numpy as np
from mpmath import fp,mp
from scipy.integrate import quad,fixed_quad, simps, dblquad
from tools.config import conf
import obslib.idis.aux

class LN:

  def __init__(self):
    #aux=conf['aux']
    #--Print out the actual values for constants

    mN     = 0.93891897
    mN2    = mN**2
    mPi2   = 0.13803
    mDel   = 1.232
    mDel2  = mDel**2
    Mbar2  = (mDel+mN)**2
    Del2   = (mDel-mN)**2
    theta  = lambda x: 0.5*(np.sign(x)+1)
    self.mN2=mN2
    self.mDel2=mDel2

    #g_PiNN = np.sqrt(14.5*4*np.pi)
    g_PiNN = np.sqrt(13.7*4*np.pi)
    gg     = g_PiNN**2/(16*np.pi**2)
    sPiN   = lambda kT2,y:(kT2+mPi2)/y+(kT2+mN2)/(1-y)
    pointN = lambda y,t,kT2: 2*gg*y/(1-y)*(-t)/(t-mPi2)**2


    g_PiND = 11.8  # in units of 1/GeV
    C      = g_PiND**2/(16*np.pi**2 * 18*mDel2)
    sPiD   = lambda kT2,y: (kT2+mPi2)/y + (kT2+mDel2)/(1-y) 
    pointD = lambda y,t,kT2: -2*C*y/(1-y)*(Mbar2-mPi2)\
                             * ( (Mbar2-mPi2)*(Del2-mPi2)-(3*(Del2-mPi2)+4*mN*mDel)*(t-mPi2) )\
                             / (t-mPi2)**2

    self.L_p2pin=0.68

    self.model=conf['SPLFUNC model'] 
    if   self.model=='IMF exp':
        self.L_p2pin=1.31
        self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * np.exp(2*(mN2-sPiN(kT2,y))/self.L_p2pin**2)
    elif self.model=='cov exp':
        self.L_p2pin=0.58
        self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * np.exp(2*(t-mPi2)/self.L_p2pin**2)
    elif self.model=='Pauli-Villars':
        self.L_p2pin=
        self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * (1 - (t-mPi2)**2/(t-self.L_p2pin**2)**2)

    # factor -2 multiplying "C" from isospin:  f_{pi+ Del0} - f_{pi- Del++} = -2 f_{pi+ Del0} => -2 f_Del^(on)
    if   self.model=='IMF exp':      self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) * np.exp(2*(mN2-sPiD(kT2,y))/self.L_p2pin**2)
    elif self.model=='cov exp':      self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) * np.exp(2*(t-mPi2)/self.L_p2pin**2)
    elif self.model=='Pauli-Villars':self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) * (1 - (t-mPi2)**2/(t-self.L_p2pin**2)**2)
    # storage for F2p
    #self.storage={}

  def integrand(self,y,kT2):
    tN = -kT2/(1-y)-y**2*self.mN2/(1-y)
    tD = -kT2/(1-y)+y*self.mN2-y*self.mDel2/(1-y)
    return self.tgrandN(y,tN,kT2)+self.tgrandD(y,tD,kT2)

  def integrandN(self,y,kT2):
    tN = -kT2/(1-y)-y**2*self.mN2/(1-y)
    return self.tgrandN(y,tN,kT2)

  def integrandD(self,y,kT2):
    tD = -kT2/(1-y)+y*self.mN2-y*self.mDel2/(1-y)
    return self.tgrandD(y,tD,kT2)

  def get_f(self,y,kT2max):
    return quad(lambda kT2: self.integrand(y,kT2),0,kT2max)[0] 

  def get_fN(self,y,kT2max):
    return quad(lambda kT2: self.integrandN(y,kT2),0,kT2max)[0]

  def get_fD(self,y,kT2max):
    return quad(lambda kT2: self.integrandD(y,kT2),0,kT2max)[0]
 
  def get_F2LN(self,xpi,y,Q2,kT2max):
    F2pi=conf['pion-stfuncs'].get_FXN(xpi,Q2)
    #F2pi=self.N*xpi**(self.a0+self.a1*self.eta(Q2))*(1-xpi)
    if conf['ln mode']=='model':
      f=self.get_fN(y,kT2max)
    elif conf['ln mode']=='piecewise':
      if y not in F: F[y]=self.get_fN(y,kT2max)
      f=F[y]
    return F2pi * f

  # moments of spliting functions

  def get_fmom(self):
    return quad(lambda y: self.get_f(y,np.inf),0,1)[0]

  def get_fNmom(self):
    return quad(lambda y: self.get_fN(y,np.inf),0,1)[0]

  def get_fDmom(self):
    return quad(lambda y: self.get_fD(y,np.inf),0,1)[0]

  #-> not sure what is the purpose of the next 2 functions

  def get_multiplicities(L_p2pin):
    L_p2pin = L
    return (3.0/2.0)*quad(lambda y: self.get_f(y,np.inf),0,1)[0]

  def get_f_tdis(self,y,kT2min,kT2max):
    return quad(lambda kT2: self.integrand(y,kT2),kT2min,kT2max)[0]


if __name__=='__main__':

  conf['Q20'] = 1.0
  conf['alphaSmode']='backward'
  conf['dglap mode']='truncated'
  conf['order']='NLO'
  conf['scheme']='ZMVFS'
  conf['ln mode']='model'
  conf['SPLFUNC model']='cov exp'


  conf['aux']=qcdlib.aux.AUX()
  conf['eweak']=qcdlib.eweak.EWEAK()
  conf['alphaS']=qcdlib.alphaS.ALPHAS()
  conf['mellin-ext']=qcdlib.mellin.MELLIN(extended=True)
  conf['mellin']=qcdlib.mellin.MELLIN()
  conf['pdf-pion']=qcdlib.pdfpion0.PDF()
  conf['pion-stfuncs']=STFUNCS()
  conf['pdf']=qcdlib.pdf0.PDF()
  conf['stfuncs']=obslib.dis.stfuncs.STFUNCS()
  ln=LN()

  x=0.000224
  Q2=7.3 
  y=0.635
  xpi=x/y
  kT2max=0.4
  deltaxL=0.2
  print ln.get_F2LN(xpi,y,Q2,kT2max)
  print ln.get_R(x,xpi,deltaxL,y,Q2,kT2max)
  print ln.get_dbub(x,Q2)


