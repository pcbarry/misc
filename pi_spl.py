#!/usr/bin/env python
import numpy as np
from mpmath import fp,mp
from scipy.integrate import quad

class SPL_FUNC:

  def __init__(self,model='IMF exp'):

    mN     = 0.93891897
    mN2    = mN**2
    mPi2   = 0.13803
    theta  = lambda x: 0.5*(np.sign(x)+1)
    self.mN2=mN2

    g_PiNN = np.sqrt(13.7*4*np.pi)
    gg     = g_PiNN**2/(16*np.pi**2)
    sPiN   = lambda kT2,y:(kT2+mPi2)/y+(kT2+mN2)/(1-y)
    pointN = lambda y,t,kT2: 2*gg*y/(1-y)*(-t)/(t-mPi2)**2

    self.model=model
    if   self.model=='IMF exp':
         self.L_p2pin=1.31
         self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * np.exp(2*(mN2-sPiN(kT2,y))/self.L_p2pin**2)
    elif self.model=='cov exp':
         self.L_p2pin=0.58
         self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * np.exp(2*(t-mPi2)/self.L_p2pin**2)
    elif self.model=='Pauli-Villars':
         self.L_p2pin=0.25
         self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * (1 - (t-mPi2)**2/(t-self.L_p2pin**2)**2)

  def integrandN(self,y,kT2):
    tN = -kT2/(1-y)-y**2*self.mN2/(1-y)
    return self.tgrandN(y,tN,kT2)

  def get_fN(self,y,kT2max):
    return quad(lambda kT2: self.integrandN(y,kT2),0,kT2max)[0]

if __name__=='__main__':

  spl=SPL_FUNC(model='IMF exp') #--This corresponds with "s-exponential"
  #spl=SPL_FUNC(model='cov exp') #--This corresponds with "t-exponential"
  #spl=SPL_FUNC(model='Pauli-Villars')

  y=0.635
  kT2max=0.4

  print spl.get_fN(y,kT2max)


