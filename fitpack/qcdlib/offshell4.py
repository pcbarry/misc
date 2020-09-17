#!/usr/bin/env python
from mpmath import fp
import numpy as np
from scipy.integrate import fixed_quad
from scipy.special import gamma
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.mellin import MELLIN,DMELLIN
from qcdlib import pdf2
from tools.config import conf

#!!!!!!!
#Currently do not know how to set sumrules without having Q2 as an argument
#!!!!!!!


class OFFSHELL: 

  def __init__(self):
    self.aux=conf['aux']
    self.mellin=conf['mellin']
    self.dmellin=conf['dmellin']
    self.pdf=conf['pdf']

    self.set_params()
    self.setup()
    
  def set_params(self):
    """
    offshell(x) = norm*(x-x0)(x-x1)(1+x0-x)
    """
    params={}
    params['proton']  =np.array([8.0,0.5,0.0])
    params['neutron'] =np.array([1.0,0.5,0.0])
    self.params=params

  def set_sumrules(self,Q2):
    #print len(self.storage.keys())
    y,w = np.polynomial.legendre.leggauss(20)
    z = 0.5*(y + 1.0)
    jac = 0.5
    if Q2 not in self.storage:
      self.storage[Q2]={}
      self.storage[Q2]['uv']={}
      self.storage[Q2]['dv']={}
      for i in range(4):
        self.storage[Q2]['uv'][i+1]=0
        self.storage[Q2]['dv'][i+1]=0
        for j in range(len(z)):
            uv = self.pdf.get_xF(z[j],Q2,'uv',evolve=True)/z[j]
            dv = self.pdf.get_xF(z[j],Q2,'dv',evolve=True)/z[j]
            int_uv = z[j]**i * uv
            int_dv = z[j]**i * dv
            self.storage[Q2]['uv'][i+1]+=w[j]*int_uv*jac
            self.storage[Q2]['dv'][i+1]+=w[j]*int_dv*jac

    v1=self.storage[Q2]['uv'][1]+self.storage[Q2]['dv'][1]
    v2=self.storage[Q2]['uv'][2]+self.storage[Q2]['dv'][2]
    v3=self.storage[Q2]['uv'][3]+self.storage[Q2]['dv'][3]
    v4=self.storage[Q2]['uv'][4]+self.storage[Q2]['dv'][4]

    for flav in self.params:
        x0=self.params[flav][1]
        x1=(-v4+(2*x0+1)*v3-(x0**2+x0)*v2)/(-v3+(2*x0+1)*v2-(x0**2+x0)*v1)
        self.params[flav][2]=x1
        
  def get_moments(self,flav,M=None):
    """
    if M==None: then parametrization is to be use to compute moments along mellin contour
    else the Mth moment is returned
    """
    if np.array(M)[0]==None: M=self.dmellin.M
    norm,x0,x1=self.params[flav]
    return norm*((-1/(M+3))+((2*x0+x1+1)/(M+2))-((x0**2+x0+2*x0*x1+x1)/(M+1))+((x0**2+x0)*x1)/M)

  def setup(self):
    self.get_BC()
    self.storage={}

  def get_state(self):
    return (self.off)

  def set_state(self,state):
    self.off = state
    self.storage = {}

  def get_BC(self):
    self.off={}
    self.off['proton'] =self.get_moments('proton',self.dmellin.M)
    self.off['neutron']=self.get_moments('neutron',self.dmellin.M)

  def beta(self,a,b):
      return gamma(a)*gamma(b)/gamma(a+b)

  def test(self,x,Q2,nucleon):
    self.pdf.setup()
    self.pdf.evolve(Q2)
    self.set_sumrules(Q2)
    self.setup()
    norm,x0,x1=self.params[nucleon]
    func = lambda y: norm*(y-x0)*(y-x1)*(1+x0-y)
    print 'At Q2 = %5.4f:'%Q2
    print 'x1 = %s'%x1

    #--test that sum rule is satisfied with parameterized pdf form

    Nu,au,bu,cu,du=self.pdf.params['uv1']
    Nd,ad,bd,cd,dd=self.pdf.params['dv1']
    normu = self.beta(2+au,bu+1)+cu*self.beta(2+au+0.5,bu+1)+du*self.beta(2+au+1,bu+1)
    normd = self.beta(2+ad,bd+1)+cd*self.beta(2+ad+0.5,bd+1)+dd*self.beta(2+ad+1,bd+1)
    uv = lambda y: Nu*y**au*(1-y)**bu*(1 + cu*y**0.5 + du*y)/normu
    dv = lambda y: Nd*y**ad*(1-y)**bd*(1 + cd*y**0.5 + dd*y)/normd
    qv = lambda y: uv(y) + dv(y)
    integrand = lambda y: func(y) * qv(y)
    sumrule = fixed_quad(integrand,0,1.0,n=20)[0]
    print 'The sum rule is %5.10f (direct PDF form)'%sumrule

    sumrule = 0
    #--test that sum rule is satisfied with inverse Mellin transform of pdf
    y,w = np.polynomial.legendre.leggauss(20)
    z = 0.5*(y + 1.0)
    jac = 0.5
    for j in range(len(z)):
        _qv = (self.pdf.get_xF(z[j],Q2,'uv')+self.pdf.get_xF(z[j],Q2,'dv'))/z[j]
        f  = func(z[j])
        sumrule += w[j]*f*_qv*jac

    print 'The sum rule is %5.10f (Mellin transform of PDF)'%sumrule

    #--test that function matches inverse mellin

    print 'xF = %5.6f (direct form)'%(x*func(x))
    print 'xF = %5.6f (Mellin transform)'%(x*self.mellin.invert(x,self.off[nucleon]))

  def get_xF(self,x,nucleon):
    return x*self.mellin.invert(x,self.off[nucleon])

if __name__=='__main__':

    conf['order']='NLO'
    conf['Q20'] = 1.27**2
    conf['aux']=AUX()
    conf['dglap mode']='truncated'
    conf['mellin']=MELLIN(npts=8)
    conf['dmellin']=DMELLIN(nptsN=8,nptsM=8)
    conf['alphaS']=ALPHAS()
    conf['pdf']=pdf2.PDF()

    offshell=OFFSHELL()

    offshell.test(0.8,1.27**2,'neutron')

    offshell.test(0.8,20,'neutron')





