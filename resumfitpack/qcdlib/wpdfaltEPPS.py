#!/usr/bin/env python
import sys
from scipy.special import gamma
import numpy as np
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.dglap import DGLAP
from qcdlib.kernels import KERNELS
from qcdlib.mellin import MELLIN
from tools.config import conf
from obslib.dy.dyresmellin import DYMELLIN
from scipy.integrate import quad
import copy
from timeit import default_timer as timer

class PDF: 

  def __init__(self,mellin=None):

    """
    This PDF should be used for resummation purposes.
    Assumes pion BEAM incident on target.
    """

    self.spl='upol'
    self.Q20=conf['Q20']
    self.mc2=conf['aux'].mc2
    self.mb2=conf['aux'].mb2

    mell=copy.deepcopy(conf['mellin'])
    mell1=copy.deepcopy(conf['dymell1'])
    mell2=copy.deepcopy(conf['dymell2'])
    mell3=copy.deepcopy(conf['dymell3'])

    self.Nsize=mell.N.size #--from M-independent contour
    self.resN1size=mell1.Nmmin.size #--one contour for -M, full for +M
    self.resN2size=mell2.Nm.size #--improvements based on M value should be made
    self.resN3size=mell3.Nmmin.size #--one contour for +M, full for -M

    self.M=conf['dymell1'].M
    self.Nfullsize=len(conf['dymell1'].N)

    mell.N=np.append(mell.N,mell1.Nmmin)
    mell.N=np.append(mell.N,mell2.Nm)
    mell.N=np.append(mell.N,mell3.Nmmin)

    print 'initializing the kernel'
    self.kernel=KERNELS(mell,self.spl)
    print 'initializing dglap'
    self.dglap=DGLAP(mell,conf['alphaS'],self.kernel,conf['dglap mode'],conf['order'])
    self.mellin=mell

    self.set_params()
    self.setup()

  def set_params(self):
    """
    THESE VALUES RESULTED FROM A FIT TO EPPS16
    PROTON
    f(x) = x**(a) * (1-x)**b * (cp0+cp1*x**0.5+cp2*x+cp3*x**1.5+cp4*x**2)
    array in order of parameters as shown above from left to right

    NUCLEUS
    R(x) = A + a1 * (x-xa)**2 + b1*x**(10*xa) + b2*x**(20*xa) + b3*x**(20*xa)\
             + (c1-c2*x)*(1-x)**(c3)
    f_p^W(x) = R(x) * f_p(x)
    up:   xf(x) = x/A*(Z*u_p^W(x)+(A-Z)*d_p^W(x))
    down: xf(x) = x/A*(Z*d_p^W(x)+(A-Z)*u_p^W(x))
    
    arrays in order of: xa,A,a1,b1,b2,b3,c1,c2,c3,a,b,cp0,cp1,cp2,cp3,cp4
    len(array)=15
    """
    params={}


    #params['g'] =np.array([6.37591335e-3,8.85313100e+0,3.56723352e+0,\
    #                        -5.49813727e+0,-7.35515980e+0,3.09812586e+0,\
    #                        3.11559299e+0,3.23288827e+1,1.74397844e+1,\
    #                        3.49919010e-1,9.66714573e+0,4.45945823e-1,\
    #                        7.04583740e+0,7.93333767e+0,1.40215881e+2,\
    #                        8.20381120e+2])
    params['g']=np.array([1.39748138e-01,5.49396534e+00,7.81331000e+0,8.57318083e-01,\
                            4.50484910e+01,1.67540378e+01,-1.14539165e+01,-3.12305023e+02,\
                            3.32640719e+01,-8.03415819e-02,6.51039283e+00,1.35199641e-02,\
                            1.35780306e+00,2.25696769e+00,-1.00531369e+01,7.49157186e+00])
    #params['uv']=np.array([0.67935427,1.00510332,7.57161783,-2.19612347,\
    #                        2.96349237,0.82632013,-3.00143803,7.29135935,\
    #                        17.35317963,0.65100653,3.37957101,1.89340062,\
    #                        -4.91606,8.77534342,0.1714838,5.04766934])
    params['uv']=np.array([1.39652505,0.75468464,-0.44925927,3.63130059,-7.36708953,24.82834897,\
                            0.10372526,-3.10184221,3.7592546,-0.29698974,3.77140151,1.28680689,\
                           -1.69814291,5.6984256,6.1093531,6.61335472])
    #params['dv']=np.array([8.99874693e-1,2.19719246e+0,-1.93603030e+0,\
    #                        5.41624014e+0,-3.68765001e+1,1.45779661e+2,\
    #                        4.98662130e-3,2.64666428e+0,3.56420363e-1,\
    #                        1.77503456e-1,3.80249682e+0,1.19666991e-2,\
    #                        2.55960155e+0,2.01313875e+0,1.66858304e-1,\
    #                        -4.69980749e+0])
    params['dv']=np.array([4.35562972e+00,-3.64734002e+00,3.67470251e-01,3.30612131e+01,\
                           -3.97654314e+02,9.74215858e+02,1.89410135e+00,-9.98227868e+00,\
                            3.60228170e+00,6.06072380e-01,2.91109547e+00,2.45085976e-01,\
                            2.39830993e-01,2.89919345e-02,8.56733354e-01,-1.09874426e+00])
    #params['ub']=np.array([4.54178018e-1,2.08046736e-1,2.98020415e+0,\
    #                        -3.11552642e+1,3.23664361e+2,-9.37369617e+2,\
    #                        -7.93897627e-1,-4.00432286e+0,9.59869678e-1,\
    #                        -1.87735773e-1,7.01821775e+0,2.25054027e+0,\
    #                        -1.20861141e+1,4.28919939e+1,-8.60313042e+1,\
    #                        6.67542167e+1])
    params['ub']=np.array([0.65713903,0.35948647,1.96038902,-2.4933741,14.38947031,\
                           -30.16467897,1.75106812,9.7922718,9.9481934,0.68585452,\
                             5.89100752,0.75145841,0.87120026,-1.66966697,-4.06491063,5.40509331])
    #params['db']=np.array([9.97044213e-1,2.87793415e+0,-2.42806990e+0,\
    #                        -5.03181881e+1,-1.99537116e+2,-7.08364532e+3,\
    #                        -1.79129725e-1,6.14635174e+0,6.67702437e-1,\
    #                        -3.33163861e-1,1.12096789e+1,5.37062437e-2,\
    #                        7.05459153e-2,1.10368424e+1,-3.95332073e+1,\
    #                        1.17885909e+2])
    params['db']=np.array([1.76065943e+00,3.24778824e-01,2.13690382e-01,1.56956857e+05,\
                            1.83971876e+07,1.81290877e+09,1.91795593e-01,1.56101603e-01,\
                           -7.06665116e+00,2.90891023e-01,1.31619296e+01,3.24214268e-01,\
                            1.02300675e+00,9.93617306e-01,4.65893356e+00,-7.16985111e+00])
    #params['s'] =np.array([8.60841465e-1,-2.35416404e-2,2.84969601e-1,\
    #                        -8.46984746e-2,3.86426949e+1,-2.44816292e+2,\
    #                        7.39546712e-1,-6.20117693e+0,1.78384159e+1,\
    #                        -1.90966762e-1,5.83833708e+0,2.42600060e-2,\
    #                        3.08044128e-2,7.63577448e-1,2.11299637e+0,\
    #                        1.74066889e+0])
    params['s']=np.array([5.36886030e-01,4.97341664e-01,-1.27649898e-01,2.68228347e+01,\
                           4.20877080e+01,2.09409713e+03,7.19791040e-01,5.49269423e+00,\
                           8.09291360e+00,8.18445569e-01,1.06329586e+01,1.45994170e+00,\
                           1.82289911e+00,2.99812775e+00,-7.99788538e+00,4.32424648e-01])

    self.params=params

  def set_sumrules(self):
    """
    NEEDS SIGNIFICANT WORK HERE
    """
    #--valence
    self.params['ubv'][0]=1    
    self.params['ubv'][0]=1/self.get_moments('ubv',1)

    self.params['dv'][0]=1    
    self.params['dv'][0]=1/self.get_moments('dv',1)

    #--msr
    up=2*self.get_moments('u',2) + self.get_moments('ubv',2)
    dp=self.get_moments('dv',2) + 2*self.get_moments('db',2)
    sp=self.get_moments('s',2)+self.get_moments('sb',2)
    self.params['g'][0]=1    
    self.params['g'][0]=(1-up-dp-sp)/self.get_moments('g',2)
    g=self.get_moments('g',2)
    msr=g+up+dp+sp

  def set_moms(self):

    ub=self.get_moments('ub')
    db=self.get_moments('db')
    s =self.get_moments('s')
    up=self.get_moments('uv')+2*ub
    dp=self.get_moments('dv')+2*db

    moms={}

    moms['g']=self.get_moments('g')
    moms['up']=conf['Z']/conf['A']*(up)+(conf['A']-conf['Z'])/conf['A']*(dp)
    moms['dp']=conf['Z']/conf['A']*(dp)+(conf['A']-conf['Z'])/conf['A']*(up)
    moms['sp']=2*s
    moms['um']=conf['Z']/conf['A']*(self.get_moments('uv'))+(conf['A']-conf['Z'])/conf['A']*(self.get_moments('dv'))
    moms['dm']=conf['Z']/conf['A']*(self.get_moments('dv'))+(conf['A']-conf['Z'])/conf['A']*(self.get_moments('uv'))
    moms['sm']=s-s
    self.moms0=moms
    self.get_BC(moms)

  def setup(self):
  
    #self.set_sumrules()
    self.set_moms()

    # we will store all Q2 values that has been precalc
    self.storage={}
    self.res_storage={}

  def beta(self,a,b):
    return gamma(a)*gamma(b)/gamma(a+b)

  def sum_cm(self,a,b,cp0,cp1,cp2,cp3,cp4):
    return cp0*self.beta(a+0.0/2,b)\
          +cp1*self.beta(a+1.0/2,b)\
          +cp2*self.beta(a+2.0/2,b)\
          +cp3*self.beta(a+3.0/2,b)\
          +cp4*self.beta(a+4.0/2,b)

  def get_moments(self,flav,N=None):
    """
    if N==None: then parametrization is to be use to compute moments along mellin contour
    else the Nth moment is returned
    """
    if N==None: N=self.mellin.N
    xa,A,a1,b1,b2,b3,c1,c2,c3,alpha,beta,cp0,cp1,cp2,cp3,cp4=self.params[flav]
    termA=np.array([A*self.sum_cm(n+alpha-1,beta+1,cp0,cp1,cp2,cp3,cp4) for n in N])
    terma1=np.array([a1*(self.sum_cm(n+alpha-1+2,beta+1,cp0,cp1,cp2,cp3,cp4)\
              -2*xa*self.sum_cm(n+alpha-1+1,beta+1,cp0,cp1,cp2,cp3,cp4)\
              +xa**2*self.sum_cm(n+alpha-1,beta+1,cp0,cp1,cp2,cp3,cp4)) for n in N])
    termb1=np.array([b1*self.sum_cm(n+alpha-1+10*xa,beta+1,cp0,cp1,cp2,cp3,cp4) for n in N])
    termb2=np.array([b2*self.sum_cm(n+alpha-1+20*xa,beta+1,cp0,cp1,cp2,cp3,cp4) for n in N])
    termb3=np.array([b3*self.sum_cm(n+alpha-1+30*xa,beta+1,cp0,cp1,cp2,cp3,cp4) for n in N])
    termc1=np.array([c1*self.sum_cm(n+alpha-1,beta+1+c3,cp0,cp1,cp2,cp3,cp4) for n in N])
    termc2=np.array([-c2*self.sum_cm(n+alpha-1+1,beta+1+c3,cp0,cp1,cp2,cp3,cp4) for n in N])
    mom=termA+terma1+termb1+termb2+termb3+termc1+termc2
    
    return mom

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
      self.res_storage={}

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
      self.storage[Q2]={}
      self.res_storage[(-1,Q2)]={}
      self.res_storage[(-2,Q2)]={}
      self.res_storage[(-3,Q2)]={}

      if self.mb2<Q2:
        moms=self.dglap.evolve(self.BC5,self.mb2,Q2,5)
      elif self.mc2<=Q2 and Q2<=self.mb2:
        moms=self.dglap.evolve(self.BC4,self.mc2,Q2,4)
      elif Q2<self.mc2:
        moms=self.dglap.evolve(self.BC3,self.Q20,Q2,3)

      i0=self.Nsize
      i1=i0+self.resN1size
      i2=i1+self.resN2size
      flav=['g','u','ub','d','db','s','sb','c','cb','b','bb']
      for f in flav:
        self.storage[Q2][f]=moms[f][:i0]
        self.res_storage[(-1,Q2)][f]=moms[f][i0:i1]
        self.res_storage[(-2,Q2)][f]=moms[f][i1:i2]
        self.res_storage[(-3,Q2)][f]=moms[f][i2:]

      #--Here, I repack the full contours in the moments
      M=self.M
      shape1=self.Nfullsize/len(M)
      shape=(shape1,self.resN1size/shape1)
      for f in flav:
        res_1=self.res_storage[(-1,Q2)][f].reshape(shape)
        repeated_1=np.einsum('i,j->ij',res_1.T[0],np.ones(len(M)/2-1))
        self.res_storage[(-1,Q2)][f]=np.append(repeated_1,res_1,axis=1).flatten()
        res_3=self.res_storage[(-3,Q2)][f].reshape(shape)
        repeated_3=np.einsum('i,j->ij',res_3.T[-1],np.ones(len(M)/2-1))
        self.res_storage[(-3,Q2)][f]=np.append(res_3,repeated_3,axis=1).flatten()

  def get_xF(self,x,Q2,flav,evolve=True):
    if evolve: self.evolve(Q2)
    return x*conf['mellin'].invert(x,self.storage[Q2][flav])
    #return x*self.mellin.invert(x,self.storage[Q2][flav])

  def get_xF0(self,x,flav,evolve=True):
    if flav=='um': mom=self.moms0['um']
    elif flav=='dm': mom=self.moms0['dm']
    elif flav=='sm': mom=self.moms0['sm']
    return x*self.mellin.invert(x,mom)

  def get_pdfs(self,x,Q2,evolve=True):
    xf=[]
    for f in ['g','u','ub','d','db','s','sb','c','cb','b','bb']:
      xf.append(self.get_xF(x,Q2,f))
    return xf/x

if __name__=='__main__':

  conf['order']='NLO'
  conf['Q20'] = 1.3**2
  conf['aux']=AUX()
  conf['mellin']=MELLIN()
  conf['dglap mode']='truncated'
  conf['alphaS']=ALPHAS()
  conf['Z']=74.
  conf['A']=184.
  #pdf=PDF(mellin=DYMELLIN(1),sign=1)
  pdf=PDF()
  #pdf.evolve(10.)
  print pdf.get_xF(0.5,10.,'b')


