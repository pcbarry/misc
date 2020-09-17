#!/bin/env python
import sys
from scipy.special import gamma
import numpy as np
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.dglap import DGLAP
from qcdlib.kernels import KERNELS
from qcdlib.mellin import MELLIN
from tools.config import conf
from scipy.integrate import quad

class PDF: 

  def __init__(self,mellin=None):

    self.spl='upol'
    self.Q20=conf['Q20']
    self.mc2=conf['aux'].mc2
    self.mb2=conf['aux'].mb2

    if mellin==None:
      self.kernel=KERNELS(conf['mellin-pion'],self.spl)
      self.dglap=DGLAP(conf['mellin-pion'],conf['alphaS'],self.kernel,conf['dglap mode'],conf['order'])
      self.mellin=conf['mellin-pion']
    else:
      self.kernel=KERNELS(mellin,self.spl)
      self.glap=DGLAP(mellin,conf['alphaS'],self.kernel,conf['dglap mode'],conf['order'])
      self.mellin=mellin
    #if 'mode' in conf: mode=conf['mode']
    #else: mode='truncated'
    self.set_params()
    self.setup()

  def set_params(self):
    """
    f(x) = norm * x**a * (1-x)**b * (1+c*x**0.5+d*x)
    """
    params={}

    params['g1'] =np.array([0,-0.5,3,0,0])
    params['g2'] =np.array([0,-0.5,3,0,0])
    params['g3'] =np.array([0,-0.5,3,0,0])
    params['ubv1']=np.array([1, 0.5,3,0,0])
    params['ubv2']=np.array([0, 0.5,3,0,0])
    params['ubv3']=np.array([0, 0.5,3,0,0])
    params['dv1']=np.array([1, 0.5,3,0,0])
    params['dv2']=np.array([0, 0.5,3,0,0])
    params['dv3']=np.array([0, 0.5,3,0,0])
    params['u1']=np.array([1,-0.5,3,0,0])
    params['u2']=np.array([0,-0.5,3,0,0])
    params['u3']=np.array([0,-0.5,3,0,0])
    params['db1']=np.array([1,-0.5,3,0,0])
    params['db2']=np.array([0,-0.5,3,0,0])
    params['db3']=np.array([0,-0.5,3,0,0])
    params['s1'] =np.array([0,-0.5,3,0,0])
    params['s2'] =np.array([0,-0.5,3,0,0])
    params['s3'] =np.array([0,-0.5,3,0,0])
    params['sb1']=np.array([0,-0.5,3,0,0])
    params['sb2']=np.array([0,-0.5,3,0,0])
    params['sb3']=np.array([0,-0.5,3,0,0])

    self.params=params

  def set_sumrules(self):


    #--valence
    self.params['ubv1'][0]=1    
    self.params['ubv1'][0]=(1-self.get_moments('ubv2',1)-self.get_moments('ubv3',1))/self.get_moments('ubv1',1)

    self.params['dv1'][0]=1    
    self.params['dv1'][0]=(1-self.get_moments('dv2',1)-self.get_moments('dv3',1))/self.get_moments('dv1',1)

    #--msr
    u23   =                            self.get_moments('u2',2)+self.get_moments('u3',2)
    ubv   = self.get_moments('ubv1',2)+self.get_moments('ubv2',2)+self.get_moments('ubv3',2)
    up23  = 2*u23 + ubv

    db23 =                           self.get_moments('db2',2)+self.get_moments('db3',2)
    dv   = self.get_moments('dv1',2)+self.get_moments('dv2',2)+self.get_moments('dv3',2)
    dp23 = dv + 2*db23

    s23  =                           self.get_moments('s2',2)+self.get_moments('s3',2)
    sb23 =                           self.get_moments('sb2',2)+self.get_moments('sb3',2)
    sp23 = s23+sb23   

    g    = self.get_moments('g1',2)+self.get_moments('g2',2)+self.get_moments('g3',2)

    self.params['u1'][0]=6    
    self.params['u1'][0]=(1-up23-dp23-sp23-g)/self.get_moments('u1',2)

    self.params['db1'][0]=self.params['u1'][0]
    self.params['s1'][0] =self.params['u1'][0]
    self.params['sb1'][0]=self.params['u1'][0]

    u  = self.get_moments('u1',2)+u23
    db = self.get_moments('db1',2)+db23
    s  = self.get_moments('s1',2)+s23
    sb = self.get_moments('sb1',2)+sb23

    up = 2*u + ubv
    dp = dv + 2*db
    sp = s + sb

    msr=g+up+dp+sp

    #--share
    self.sr={}
    self.sr['msr']=msr
    self.sr['ubv(1)'] = self.get_moments('ubv1',1)+self.get_moments('ubv2',1)+self.get_moments('ubv3',1)
    self.sr['dv(1)']  = self.get_moments('dv1',1)+self.get_moments('dv2',1)+self.get_moments('dv3',1)
    self.sr['g(2)']   = g
    self.sr['sea(2)'] = 6*u 
    self.sr['valence(2)']=2*ubv

  def set_moms(self):

    moms={}

    moms['g']=self.get_moments('g1')+self.get_moments('g2')+self.get_moments('g3')

    u   = self.get_moments('u1')+self.get_moments('u2')+self.get_moments('u3')
    ubv = self.get_moments('ubv1')+self.get_moments('ubv2')+self.get_moments('ubv3')
    ub  = ubv+u
    moms['up']=u + ub

    db = self.get_moments('db1')+self.get_moments('db2')+self.get_moments('db3')
    dv = self.get_moments('dv1')+self.get_moments('dv2')+self.get_moments('dv3')
    d  = dv+db
    moms['dp']=d+db

    s  = self.get_moments('s1') +self.get_moments('s2')+self.get_moments('s3')
    sb = self.get_moments('sb1')+self.get_moments('sb2')+self.get_moments('sb3')
    moms['sp']=s+sb   

    moms['um']=u-ub
    moms['dm']=d-db
    moms['sm']=s-sb
    self.moms0=moms
    self.get_BC(moms)

  def setup(self):
  
    self.set_sumrules()
    self.set_moms()

    # we will store all Q2 values that has been precalc
    self.storage={}

  def beta(self,a,b):
    return gamma(a)*gamma(b)/gamma(a+b)

  def get_moments(self,flav,N=None):
    """
    if N==None: then parametrization is to be use to compute moments along mellin contour
    else the Nth moment is returned
    """
    if N==None: N=conf['mellin-pion'].N
    M,a,b,c,d=self.params[flav]
    #if 'v' in flav:
    if 'resum' in conf['pdf-pion parametrization'] and 'v' in flav:
        norm=self.beta(2+a,b+1)+c*self.beta(2+a+d,b+1)
        mom=self.beta(N+a,b+1)+c*self.beta(N+a+d,b+1)
    else:
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

    ###############################################
    # BC for Nf=3
    g   = moms['g']
    up  = moms['up']
    um  = moms['um']
    dp  = moms['dp']
    dm  = moms['dm']
    sp  = moms['sp']
    sm  = moms['sm']
    cp  = zero
    cm  = zero
    bp  = zero
    bm  = zero
    self.BC3=self._get_BC(g,up,um,dp,dm,sp,sm,zero,zero,zero,zero,zero,zero)

    ###############################################
    # BC for Nf=4
    BC4=self.dglap.evolve(self.BC3,self.Q20,self.mc2,3)
    g =BC4['g']
    up=BC4['up']
    dp=BC4['dp']
    sp=BC4['sp']
    cp=BC4['cp']
    bp=BC4['bp']
    tp=BC4['tp']
    um=BC4['um']
    dm=BC4['dm']
    sm=BC4['sm']
    cm=BC4['cm']
    bm=BC4['bm']
    tm=BC4['tm']
    self.BC4=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,tp,tm)

    ###############################################
    # BC for Nf=5
    BC5=self.dglap.evolve(self.BC4,self.mc2,self.mb2,4)
    g =BC5['g']
    up=BC5['up']
    dp=BC5['dp']
    sp=BC5['sp']
    cp=BC5['cp']
    bp=BC5['bp']
    tp=BC5['tp']
    um=BC5['um']
    dm=BC5['dm']
    sm=BC5['sm']
    cm=BC5['cm']
    bm=BC5['bm']
    tm=BC5['tm']
    self.BC5=self._get_BC(g,up,um,dp,dm,sp,sm,cp,cm,bp,bm,tp,tm)

  def evolve(self,Q2):

    if Q2 not in self.storage:
      if self.mb2<Q2: 
        self.storage[Q2]=self.dglap.evolve(self.BC5,self.mb2,Q2,5)
      elif self.mc2<=Q2 and Q2<=self.mb2: 
        self.storage[Q2]=self.dglap.evolve(self.BC4,self.mc2,Q2,4)
      elif Q2<self.mc2: 
        self.storage[Q2]=self.dglap.evolve(self.BC3,self.Q20,Q2,3)
  
  def get_xF(self,x,Q2,flav,evolve=True):
    if evolve: self.evolve(Q2)
    #return x*conf['mellin-pion'].invert(x,self.storage[Q2][flav])
    return x*self.mellin.invert(x,self.storage[Q2][flav])

  def get_xF0(self,x,flav,evolve=True):
    if flav=='um': mom=self.moms0['um']
    elif flav=='dm': mom=self.moms0['dm']
    elif flav=='sm': mom=self.moms0['sm']
    #return x*conf['mellin-pion'].invert(x,mom)
    return x*self.mellin.invert(x,mom)

  def get_second_moment(self,flav,Q2):
    self.evolve(Q2)
    if flav=='g': 
      f=lambda x: self.get_xF(x,Q2,'g',evolve=False)
    if flav=='sea': 
      f=lambda x: 2*self.get_xF(x,Q2,'u',evolve=False)\
                 +2*self.get_xF(x,Q2,'db',evolve=False)\
                 +self.get_xF(x,Q2,'s',evolve=False)\
                 +self.get_xF(x,Q2,'sb',evolve=False)
    if flav=='valence':
      f=lambda x: self.get_xF(x,Q2,'ub',evolve=False)-self.get_xF(x,Q2,'u',evolve=False)\
                 +self.get_xF(x,Q2,'d',evolve=False) -self.get_xF(x,Q2,'db',evolve=False)
    return quad(f,0,1)[0]

  def get_pdfs(self,x,Q2,evolve=True):
    xf=[]
    for f in ['g','u','ub','d','db','s','sb','c','cb','b','bb']:
      xf.append(self.get_xF(x,Q2,f))
    return xf/x

if __name__=='__main__':

  conf['order']='NLO'
  conf['Q20'] = 1.0
  conf['aux']=AUX()
  conf['mellin']=MELLIN(npts=8)
  conf['alphaS']=ALPHAS()
  pdf=PDF('upol')
  ppdf=PDF('pol')
  print pdf.get_xF(0.5,10.,'u')
  print ppdf.get_xF(0.5,10.,'u')


