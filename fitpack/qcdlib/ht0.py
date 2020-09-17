#!/usr/bin/env python
from mpmath import fp
import numpy as np
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.mellin import MELLIN
from tools.config import conf

class T3PPDF:
  
    def __init__(self):
  
        self.aux=conf['aux']
        self.Q20=conf['Q20']
        self.mc2=self.aux.mc2
        self.mb2=self.aux.mb2
        self.mellin=conf['mellin']
        conf['alphaS']=conf['alphaS']
        self.setup_kernel()
        self.set_params()
        self.setup()
  
    def setup_kernel(self):
        # anomalous dimension
        N=self.mellin.N
        psi=np.array([fp.psi(0,complex(n.real,n.imag)) for n in N])
        self.gamma=2*3*(psi+fp.euler-1/4.0+1/N/2)
  
    def set_params(self):
        """
        PDF(x) = norm * x**a * (1-x)**b * (1+c*x**0.5+d*x)
        """
        params={}
        params['u']=np.array([1, 0.5,3,0,0])
        params['d']=np.array([1, 0.5,3,0,0])
        self.params=params
  
    def setup(self):
        self.get_BC()
        self.storage={}
  
    def get_moments(self,flav,N=None):
        """
        if N==None: then parametrization is to be use to compute moments along mellin contour
        else the Nth moment is returned
        """
        if N==None: N=self.mellin.N
        norm,a,b,c,d=self.params[flav]
        return norm*(self.aux.beta(N+a,b+1)+c*self.aux.beta(N+a+0.5,b+1)+d*self.aux.beta(N+a+1,b+1))
  
    def get_BC(self):
  
        self.BC3={}
        self.BC3['u']=self.get_moments('u')
        self.BC3['d']=self.get_moments('d')
        self.BC4=self._evolve(self.BC3,self.Q20,self.mc2,3)
        self.BC5=self._evolve(self.BC4,self.mc2,self.mb2,4)
  
    def largeXevolve(self,Q20,Q2,D0,Nf):
        a=conf['alphaS'].get_a(Q2)
        a0=conf['alphaS'].get_a(Q20)
        b=11.0-2/3.*Nf
        return (a/a0)**(self.gamma/b)*D0
  
    def _evolve(self,D0,Q2a,Q2b,Nf):
        out={}
        out['u']=self.largeXevolve(Q2a,Q2b,D0['u'],Nf)
        out['d']=self.largeXevolve(Q2a,Q2b,D0['d'],Nf)
        return out
  
    def evolve(self,Q2):
  
        if Q2 not in self.storage:
            if self.mb2<Q2: 
                self.storage[Q2]=self._evolve(self.BC5,self.mb2,Q2,5)
            elif self.mc2<=Q2 and Q2<=self.mb2: 
                self.storage[Q2]=self._evolve(self.BC4,self.mc2,Q2,4)
            elif Q2<self.mc2: 
                self.storage[Q2]=self._evolve(self.BC3,self.Q20,Q2,3)
  
    def get_xF(self,x,Q2,flav,evolve=True):
        if evolve: self.evolve(Q2)
        return x*self.mellin.invert(x,self.storage[Q2][flav])
  
class T4: 

  def __init__(self):
    self.aux=conf['aux']
    self.mellin=conf['mellin']
    self.set_params()
    self.setup()

  def set_params(self):
    params={}
    params['proton']  =np.array([1,0.5,3,0,0])
    params['neutron'] =np.array([1,0.5,3,0,0])
    self.params=params

  def get_moments(self,flav,N=None):
    """
    if N==None: then parametrization is to be use to compute moments along mellin contour
    else the Nth moment is returned
    """
    if np.array(N)[0]==None: N=self.mellin.N
    M,a,b,c,d=self.params[flav]
    norm = self.aux.beta(2+a,b+1)+c*self.aux.beta(2+a+0.5,b+1)+d*self.aux.beta(2+a+1,b+1)
    mom = self.aux.beta(N+a,b+1)+c*self.aux.beta(N+a+0.5,b+1)+d*self.aux.beta(N+a+1,b+1)
    return M*mom/norm

  def setup(self):
    self.get_BC()

  def get_state(self):
    return (self.H)

  def set_state(self,state):
    self.H = state
    self.storage = {}

  def get_BC(self):
    self.H={}
    self.H['proton'] =self.get_moments('proton',self.mellin.N)
    self.H['neutron']=self.get_moments('neutron',self.mellin.N)

  def get_xF(self,x,Q2,nucleon):
    return x*self.mellin.invert(x,self.H[nucleon])/Q2

if __name__=='__main__':

    conf['order']='NLO'
    conf['Q20'] = 1.0
    conf['aux']=AUX()
    conf['mellin']=MELLIN(npts=8)
    conf['alphaS']=ALPHAS()

    t3ppdf=T3PPDF()
    t4=T4()

    print t3ppdf.get_xF(0.5,10.,'u')
    print t4.get_xF(0.5,10.,'proton')






