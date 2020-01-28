#!/usr/bin/env python
import numpy as np
from mpmath import fp,mp
from scipy.integrate import quad,fixed_quad, simps, dblquad
from tools.config import conf
import obslib.idis.aux

class PI_STFUNCS(obslib.idis.aux.AUX):
  '''
  DIS structure functions for pion. It is a simplified 
  copy from the nucleon DIS stfuncs. The main change is 
  the call to pion pdfs.
  '''

  def __init__(self):

    self.mellin=conf['mellin-pion']
    if conf['order']=='LO': self.order=0
    if conf['order']=='NLO': self.order=1
    self.setup()

  # twist level unpolarized structure functions in Mellin space

  def get_T2CFX(self,stf,nucleon,Q2,evolve=True):  
    """
    CF(2,L,3) = F2/x,FL/x,F3  
    """ 
    if evolve: conf['pdf-pion'].evolve(Q2)

    g =conf['pdf-pion'].storage[Q2]['g']     # gluon[N] were N mellin contour
    Nf=conf['alphaS'].get_Nf(Q2) # num of active flav
    a=conf['alphaS'].get_a(Q2)   # a=alphaS/4pi

    if stf=='F2':
      CQ = self.C2Q[0] + a*self.order*self.C2Q[1]  # hard kernels 
      CG = self.C2G[0] + a*self.order*self.C2G[1]
      q=np.copy(conf['pdf-pion'].storage[Q2]['qp'])    # (q + qb)[N]
      aX='ap'

    if nucleon=="neutron":
      qup=np.copy(q[1])
      qdn=np.copy(q[2])
      q[1]=qdn
      q[2]=qup

    FX  = np.zeros(self.mellin.N.size,dtype=complex) # initialize FX[N] array

    for i in range(1,Nf+1):
      aXval = self.get_aX(aX,i,Q2)
      FX+=aXval*(CQ*q[i] + 2*CG*g)

    return FX

  # unpolarized functions

  def get_FXN(self,x,Q2,stf='F2',twist=2,nucleon='proton',tmc=True,precalc=False,evolve=True):

    if   stf=='F2': FX= x*self.get_T2CFX('F2',nucleon,Q2,evolve)
    FX*=x**(-self.mellin.N)

    return self.mellin.invert(1,FX) # here x=1 since the factor x**(-N) is included in the moments


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


