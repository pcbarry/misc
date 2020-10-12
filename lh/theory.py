#!/usr/bin/env python
import numpy as np
from mpmath import fp,mp
from scipy.integrate import quad,fixed_quad, simps, dblquad
from tools.config import conf
import obslib.idis.aux

class STFUNCS(obslib.idis.aux.AUX):
  '''
  DIS structure functions for kaon. It is a simplified 
  copy from the nucleon DIS stfuncs. The main change is 
  the call to kaon pdfs.
  '''

  def __init__(self):

    self.mellin=conf['mellin-kaon']
    #self.mellin=conf['mellin']
    if conf['order']=='LO': self.order=0
    if conf['order']=='NLO': self.order=1
    self.setup()

  # twist level unpolarized structure functions in Mellin space

  def get_T2CFX(self,stf,nucleon,Q2,evolve=True):  
    """
    CF(2,L,3) = F2/x,FL/x,F3  
    """ 
    if evolve: conf['pdf-kaon'].evolve(Q2)
    #if evolve: conf['pdf'].evolve(Q2)

    g =np.copy(conf['pdf-kaon'].storage[Q2]['g'])     # gluon[N] were N mellin contour
    #g =conf['pdf'].storage[Q2]['g']     # gluon[N] were N mellin contour
    Nf=conf['alphaS'].get_Nf(Q2) # num of active flav
    a=conf['alphaS'].get_a(Q2)   # a=alphaS/4pi

    if stf=='F2':
      CQ = self.C2Q[0] + a*self.order*self.C2Q[1]  # hard kernels 
      CG = self.C2G[0] + a*self.order*self.C2G[1]
      q=np.copy(conf['pdf-kaon'].storage[Q2]['qp'])    # (q + qb)[N]
      #q=np.copy(conf['pdf'].storage[Q2]['qp'])    # (q + qb)[N]
      aX='ap'

    elif stf=='FL':
      CQ = a*self.order*self.CLQ[1]
      CG = a*self.order*self.CLG[1]
      q=np.copy(conf['pdf-kaon'].storage[Q2]['qp'])    # (q + qb)[N]
      aX='ap'

    elif stf=='F3':
      CQ = self.C3Q[0] + a*self.order*self.C3Q[1]
      CG = 0
      q=np.copy(conf['pdf-kaon'].storage[Q2]['qm'])    # (q - qb)[N]
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

#L_p2kan=0.68  # UV regulator

#initial=True  # use this flag to cp the params to conf for the initial eval
              # see parman.set_p_to_ka_n_params to see how it works

#F={}          # dictionary to store f(y) values 

class LH:

  def __init__(self):
    aux=conf['aux']
    """
    Formalism from PRD 94, 094035 (2016)
    """
    self.mN   = aux.M
    self.mN2  = aux.M2
    self.mK   = aux.Mk
    self.mK2  = self.mK**2
    self.ML   = 1.115683 #--Lambda mass particle from PDG
    self.Mbar = self.ML+self.mN
    self.DelM = self.ML-self.mN

    self.model=conf['SPLFUNC model']

    #--prefactors

    psdc = 0.093 #--pseudoscalar decay constant

    D = 0.85
    F = 0.41
    CKL2 = ((D+3*F)/2.0/3**0.5)**2 #--K^+ Lambda coefficient C_{KY}

    prefact=CKL2*self.Mbar**2/(4*np.pi*psdc)**2
    DKY = lambda xL,kT2: -(kT2+(1-xL)*self.ML**2+xL*self.mK2-(1-xL)*xL*self.mN2)/xL
    self.fKpL = lambda xL,kT2: prefact*(1-xL)*(kT2+(self.mN*(1-xL)+self.DelM)**2)/(xL**2*DKY(xL,kT2)**2)

    #--parameter
    self.L_p2KLam=1.1

  def get_regulator(self,kT2,xL):

        sKL = (kT2+self.mK2)/(1-xL)+(kT2+self.ML**2)/xL
        t = -kT2/xL - (1-xL)/xL*(self.ML**2-xL*self.mN2)
        DKY = -(kT2+(1-xL)*self.ML**2+xL*self.mK2-(1-xL)*xL*self.mN2)/xL

        if self.model=='s-dep exp':
            #reg = np.exp(2*(self.mN2-sKL)/Lambda**2)
            reg = np.exp((self.mN2-sKL)/self.L_p2KLam**2)
        elif self.model=='t-dep exp':
            reg = np.exp(2*(t-self.mK2)/self.L_p2KLam**2)
        elif self.model=='t-dep mon':
            #reg = ((Lambda**2-self.mK2)/(Lambda**2-t))**2
            reg = ((self.L_p2KLam**2-self.mK2)/(self.L_p2KLam**2-t))**2
        elif self.model=='t-dep mon':
            #reg = ((Lambda**2-self.mK2)/(Lambda**2-t))**2
            reg = ((self.L_p2KLam**2-self.mK2)/(self.L_p2KLam**2-t))**4
        elif self.model=='Regge': 
            reg = (1-xL)**(-2*t) * np.exp(2*(t-self.mK2)/self.L_p2KLam**2)
        elif self.model=='Pauli-Villars':
            reg = (1 - (t-self.mK2)**2/(t-self.L_p2KLam**2)**2)
            #mu1=0.545 #--from PRD 94, 094035 (2016)
            #reg = 1- DKY**2/(t-mu1**2)**2
        return reg

  def get_fLkT2(self,xL,kT2):
        reg=self.get_regulator(kT2,xL)
        fL = self.fKpL(xL,kT2)*reg
        return fL

  def get_fL(self,xL,kT2max):
      return quad(lambda kT2: self.get_fLkT2(xL,kT2),0,kT2max)[0]

  def get_kinfact(self,Q2,x,ye):
      eweak=conf['eweak']
      aEM=eweak.get_alpha(Q2)
      return 4.0*np.pi*aEM**2/x/Q2**2*(1.0-ye+ye**2/2.0)

  def get_dsigdxdQ2dxL(self,x,xK,y,Q2,kT2max,ye):
      F2K=conf['kaon-stfuncs'].get_FXN(xK,Q2)
      xL=1.0-np.array(y)
      f=self.get_fL(xL,kT2max)
      K=self.get_kinfact(Q2,x,ye)
      return F2K*f*K

  # get, set state

  def get_state(self):
      return (self.L_p2KLam)

  def set_state(self,state):
      self.L_p2KLam=state

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


