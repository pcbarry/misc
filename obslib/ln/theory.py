#!/usr/bin/env python
import numpy as np
from mpmath import fp,mp
from scipy.integrate import quad,fixed_quad, simps, dblquad
from tools.config import conf
import obslib.idis.aux

class STFUNCS(obslib.idis.aux.AUX):
  '''
  DIS structure functions for pion. It is a simplified 
  copy from the nucleon DIS stfuncs. The main change is 
  the call to pion pdfs.
  '''

  def __init__(self):

    self.mellin=conf['mellin-pion']
    #self.mellin=conf['mellin']
    if conf['order']=='LO': self.order=0
    if conf['order']=='NLO': self.order=1
    self.setup()

  # twist level unpolarized structure functions in Mellin space

  def get_T2CFX(self,stf,nucleon,Q2,evolve=True):  
    """
    CF(2,L,3) = F2/x,FL/x,F3  
    """ 
    if evolve: conf['pdf-pion'].evolve(Q2)
    #if evolve: conf['pdf'].evolve(Q2)

    g =np.copy(conf['pdf-pion'].storage[Q2]['g'])     # gluon[N] were N mellin contour
    #g =conf['pdf'].storage[Q2]['g']     # gluon[N] were N mellin contour
    Nf=conf['alphaS'].get_Nf(Q2) # num of active flav
    a=conf['alphaS'].get_a(Q2)   # a=alphaS/4pi

    if stf=='F2':
      CQ = self.C2Q[0] + a*self.order*self.C2Q[1]  # hard kernels 
      CG = self.C2G[0] + a*self.order*self.C2G[1]
      q=np.copy(conf['pdf-pion'].storage[Q2]['qp'])    # (q + qb)[N]
      #q=np.copy(conf['pdf'].storage[Q2]['qp'])    # (q + qb)[N]
      aX='ap'

    elif stf=='FL':
      CQ = a*self.order*self.CLQ[1]
      CG = a*self.order*self.CLG[1]
      q=np.copy(conf['pdf-pion'].storage[Q2]['qp'])    # (q + qb)[N]
      aX='ap'

    elif stf=='F3':
      CQ = self.C3Q[0] + a*self.order*self.C3Q[1]
      CG = 0
      q=np.copy(conf['pdf-pion'].storage[Q2]['qm'])    # (q - qb)[N]
      aX='am'

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
    elif stf=='FL': FX= x*self.get_T2CFX('FL',nucleon,Q2,evolve)
    elif stf=='F3': FX=   self.get_T2CFX('F3',nucleon,Q2,evolve)
    FX*=x**(-self.mellin.N)

    return self.mellin.invert(1,FX) # here x=1 since the factor x**(-N) is included in the moments

# we define global parameters here 
# for easy connection with parman.py

#L_p2pin=0.68  # UV regulator

#initial=True  # use this flag to cp the params to conf for the initial eval
              # see parman.set_p_to_pi_n_params to see how it works

#F={}          # dictionary to store f(y) values 

class LN:

  def __init__(self):
    aux=conf['aux']

    mN     = aux.M
    mN2    = aux.M2
    mPi2   = aux.Mpi2
    mDel   = aux.Mdelta
    mDel2  = aux.Mdelta2
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
    if   self.model=='kT2 cutoff':   self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) #* theta(L_p2pin**2-kT2)
    elif self.model=='IMF exp':      self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * np.exp(2*(mN2-sPiN(kT2,y))/self.L_p2pin**2)
    elif self.model=='cov exp':      self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * np.exp(2*(t-mPi2)/self.L_p2pin**2)
    elif self.model=='cov mon':      self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * ((self.L_p2pin**2-mPi2)/(self.L_p2pin**2-t))**2
    elif self.model=='Pauli-Villars':self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * (1 - (t-mPi2)**2/(t-self.L_p2pin**2)**2)
    elif self.model=='Regge cov exp':self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * y**(-2*t) * np.exp(2*(t-mPi2)/self.L_p2pin**2)
    elif self.model=='Bishari':      self.tgrandN = lambda y,t,kT2: pointN(y,t,kT2) * y**(-2*t)

    # factor -2 multiplying "C" from isospin:  f_{pi+ Del0} - f_{pi- Del++} = -2 f_{pi+ Del0} => -2 f_Del^(on)
    if   self.model=='kT2 cutoff':   self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) #* theta(L_p2pin**2-kT2)
    elif self.model=='IMF exp':      self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) * np.exp(2*(mN2-sPiD(kT2,y))/self.L_p2pin**2)
    elif self.model=='cov exp':      self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) * np.exp(2*(t-mPi2)/self.L_p2pin**2)
    elif self.model=='cov mon':      self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) * ((self.L_p2pin**2-mPi2)/(self.L_p2pin**2-t))**2
    elif self.model=='Pauli-Villars':self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) * (1 - (t-mPi2)**2/(t-self.L_p2pin**2)**2)
    elif self.model=='Regge cov exp':self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) * y**(-2*t)*np.exp(2*(t-mPi2)/self.L_p2pin**2)
    elif self.model=='Bishari':      self.tgrandD = lambda y,t,kT2: pointD(y,t,kT2) * y**(-2*t)

    self.N=1
    self.a0=1
    self.a1=1
    LQCD2=0.4**2
    Q20=1.0
    self.eta=lambda Q2: np.log(np.log(Q2/LQCD2)/np.log(Q20/LQCD2))

    # storage for F2p
    self.storage={}

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
    if self.model=='kT2 cutoff': kT2max=np.amin([self.L_p2pin**2,kT2max])
    return quad(lambda kT2: self.integrand(y,kT2),0,kT2max)[0] 

  def get_fN(self,y,kT2max):
    if self.model=='kT2 cutoff': kT2max=np.amin([self.L_p2pin**2,kT2max])
    return quad(lambda kT2: self.integrandN(y,kT2),0,kT2max)[0]

  def get_fD(self,y,kT2max):
    if self.model=='kT2 cutoff': kT2max=np.amin([self.L_p2pin**2,kT2max])
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

  def get_R(self,x,xpi,deltaxL,y,Q2,kT2max):
    F2pi=conf['pion-stfuncs'].get_FXN(xpi,Q2)
    key='%f,%f'%(x,Q2)
    if key not in self.storage:
      #self.storage[key]=conf['stfuncs'].get_FXN(x,Q2,stf='F2',tmc=False)
      self.storage[key]=conf['idis stfuncs'].get_FXN(x,Q2,stf='F2',tmc=False)
    F2p= self.storage[key]
    if conf['ln mode']=='model':
      f=self.get_fN(y,kT2max)
    elif conf['ln mode']=='piecewise':
      if y not in F: F[y]=self.get_f(y,kT2max)
      f=F[y]
    return deltaxL * F2pi / F2p * f

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

  def get_tmom(self,model,L,ymax,kT2max):
    return quad(lambda y: self.get_f(y,model,kT2max(y)),0,ymax)[0]

  # db-ub 
  
  def get_qv(self,x,Q2):
    xqv = conf['pdf-pion'].get_xF(x,Q2,'ub')-conf['pdf-pion'].get_xF(x,Q2,'u')
    return xqv/x

  def get_dbub(self,x,Q2):
    return fixed_quad(np.vectorize(lambda y: 1.0/y*self.get_f(y,np.inf)*self.get_qv(x/y,Q2)),x,1.0,n=5)[0]

  def get_dbubN(self,x,Q2):
    return fixed_quad(np.vectorize(lambda y: 1.0/y*self.get_fN(y,np.inf)*self.get_qv(x/y,Q2)),x,1.0,n=5)[0]

  def get_dbubD(self,x,Q2):
    return fixed_quad(np.vectorize(lambda y: 1.0/y*self.get_fD(y,np.inf)*self.get_qv(x/y,Q2)),x,1.0,n=5)[0]

  # get, set state

  def get_state(self):
      return (self.L_p2pin)
  
  def set_state(self,state):
      self.L_p2pin=state

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


