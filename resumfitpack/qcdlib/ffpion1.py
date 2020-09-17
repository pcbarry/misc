#!/usr/bin/env python
import numpy as np
from scipy.special import gamma
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.dglap import DGLAP
from qcdlib.kernels import KERNELS
from qcdlib.mellin import MELLIN
from tools.config import conf
      
class FF: 

  def __init__(self):

    self.spl='upol_ff'
    self.Q20=conf['Q20']
    self.mc2=conf['aux'].mc2
    self.mb2=conf['aux'].mb2
    self.mellin=conf['mellin']
    self.kernel=KERNELS(self.mellin,self.spl)
    self.dglap=DGLAP(self.mellin,conf['alphaS'],self.kernel,conf['dglap mode'],conf['order'])
    self.set_params()
    self.setup()

  def set_params(self):
     
    ## D(z) = norm * z**a * (1-z)**b * (1+c*z**0.5+d*z)
      
    params={}

    ## first shapes
    params['g1']   = np.array([0,-0.5,3,0,0])
    params['up1']   = np.array([1,-0.5,3,0,0])
    params['dp1']   = np.array([1,-0.5,3,0,0])
    params['sp1']   = np.array([1,-0.5,3,0,0])
    params['cp1']   = np.array([1,-0.5,3,0,0])
    params['bp1']   = np.array([1,-0.5,3,0,0])

    ## second shapes
    params['up2']   = np.array([0,-0.5,3,0,0])
    params['dp2']   = np.array([0,-0.5,3,0,0])

    self.params=params

  def set_sumrules(self):
    pass

  def set_moms(self):
    moms={}
    moms['g']  =self.get_moments('g1')
    moms['up1'] =self.get_moments('up1')
    moms['up2'] =self.get_moments('up2')
    moms['dp1'] =self.get_moments('dp1')
    moms['dp2'] =self.get_moments('dp2')
    moms['sp'] =self.get_moments('sp1')
    moms['cp'] =self.get_moments('cp1')
    moms['bp'] =self.get_moments('bp1')

    moms['up']=moms['up1']+moms['up2']
    moms['dp']=moms['dp1']+moms['dp2']

    self.get_BC(moms)

  def setup(self):
    self.set_sumrules()
    self.set_moms()
    self.storage={} # we will store all Q2 values that has been precalc

  def beta(self,a,b):
    return gamma(a)*gamma(b)/gamma(a+b)

  def get_moments(self,flav,N=None):
    """
    if N==None: then parametrization is to be use to compute moments along mellin contour
    else the Nth moment is returned
    """
    if N==None: N=self.mellin.N
    M,a,b,c,d=self.params[flav]
    norm=self.beta(2+a,b+1)+c*self.beta(2+a+0.5,b+1)+d*self.beta(2+a+1.0,b+1)
    mom=self.beta(N+a,b+1)+c*self.beta(N+a+0.5,b+1)+d*self.beta(N+a+1.0,b+1)
    return M*mom/norm 

  def _get_BC(self,g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,tp,tm):
    N=self.mellin.N
    
    # flav composition
    vm,vp={},{}
    vm[35]= bm + cm + dm + sm - 5*tm + um
    vm[24]= -4*bm + cm + dm + sm + um 
    vm[15]= -3*cm + dm + sm + um
    vm[8] = dm - 2*sp + 2*(-sm + sp) + um
    vm[3] = -dm + um
    vm[0] = np.zeros(N.size,dtype=complex)
    vp[0] = np.zeros(N.size,dtype=complex)
    vp[3] = -dp + up
    vp[8] = dp - 2*sp + up 
    vp[15]= -3*cp + dp + sp + up 
    vp[24]= -4*bp + cp + dp + sp + up 
    vp[35]= bp + cp + dp + sp - 5*tp + up 
    qs    = bp + cp + dp + sp + tp + up 
    qv    = bm + cm + dm + sm + tm + um 
    q     = np.zeros((2,N.size),dtype=complex)
    q[0]=np.copy(qs)
    q[1]=np.copy(g)

    BC={}
    BC['vm']=vm 
    BC['vp']=vp 
    BC['qv']=qv
    BC['q'] =q
    return BC

  def get_state(self):
     return (self.BC3,self.BC4,self.BC5)
  
  def set_state(self,state):
      self.BC3, self.BC4, self.BC5 = state[:]
      self.storage = {}

  def get_BC(self,moms):

    N=self.mellin.N
    zero=np.zeros(N.size,dtype=complex)

    ######################################
    # BC3 
    g   = moms['g']
    up  = moms['up']
    dp  = moms['dp']
    sp  = moms['sp']
    cp  = zero 
    bp  = zero
    um  = zero
    dm  = zero
    sm  = zero
    cm  = zero 
    bm  = zero
    self.BC3=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,zero,zero)

    ######################################
    # BC for Nf=4
    BC4=self.dglap.evolve(self.BC3,self.Q20,self.mc2,3)
    g =BC4['g']
    up=BC4['up']
    dp=BC4['dp']
    sp=BC4['sp']
    cp=moms['cp']
    bp=zero
    um=zero
    dm=zero
    sm=zero
    cm=zero
    bm=zero
    self.BC4=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,zero,zero)

    ######################################
    # BC for Nf=5
    BC5=self.dglap.evolve(self.BC4,self.mc2,self.mb2,4)
    g =BC5['g']
    up=BC5['up']
    dp=BC5['dp']
    sp=BC5['sp']
    cp=BC5['cp']
    bp=moms['bp']
    um=zero
    dm=zero
    sm=zero
    cm=zero
    bm=zero
    self.BC5=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,zero,zero)

  def evolve(self,Q2):

    if Q2 not in self.storage:
      if self.mb2<=Q2: 
        self.storage[Q2]=self.dglap.evolve(self.BC5,self.mb2,Q2,5)
      elif self.mc2<=Q2 and Q2<self.mb2: 
        self.storage[Q2]=self.dglap.evolve(self.BC4,self.mc2,Q2,4)
      elif Q2<self.mc2: 
        self.storage[Q2]=self.dglap.evolve(self.BC3,self.Q20,Q2,3)
  
  def get_xF(self,x,Q2,flav,evolve=True):
    if evolve: self.evolve(Q2)
    return x*self.mellin.invert(x,self.storage[Q2][flav])

  def get_xF0(self,x,flav):
    M,a,b,c,d=self.params[flav]
    norm=self.beta(2+a,b+1)+c*self.beta(2+a+0.5,b+1)+d*self.beta(2+a+1.0,b+1)
    shape=x**a*(1-x)**b*(1+c*x**0.5+d*x) 
    return x*M*shape/norm

if __name__=='__main__':

  conf['order']='NLO'
  conf['dglap mode']='truncated'
  conf['Q20'] = 1.0
  conf['aux']=AUX()
  conf['mellin']=MELLIN(npts=8)
  conf['alphaS']=ALPHAS()
  ff=FF()
  print ff.get_xF(0.5,10.,'u')





