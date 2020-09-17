#!/usr/bin/env python
from mpmath import fp
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.mellin import MELLIN,DMELLIN
from qcdlib import pdf2
from tools.config import conf

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
    params['proton']  =np.array([8.10,0.448,0.05])
    params['neutron'] =np.array([8.10,0.448,0.05])
    self.params=params

  def set_sumrules(self):
    vmom={}
    for flav in self.params:
        vmom[flav]={}
        x0=self.params[flav][1]
        for i in range(4):
            vmom[flav][i+1]=self.pdf.get_moments('uv1',N=i+1)+self.pdf.get_moments('dv1',N=i+1)

        x1=(-vmom[flav][4]+(2*x0+1)*vmom[flav][3]-(x0**2+x0)*vmom[flav][2])/(-vmom[flav][3]+(2*x0+1)*vmom[flav][2]-(x0**2+x0)*vmom[flav][1])
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
    self.set_sumrules()
    self.get_BC()

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
    norm,x0,x1=self.params[nucleon]
    func = lambda y: norm*(y-x0)*(y-x1)*(1+x0-y)
    self.pdf.setup()
    self.pdf.evolve(Q2)
    #print x1
    print 'At Q2 = %5.4f:'%Q2

    #--test that sum rule is satisfied with parameterized pdf form

    Nu,au,bu,cu,du=self.pdf.params['uv1']
    Nd,ad,bd,cd,dd=self.pdf.params['dv1']
    normu = self.beta(2+au,bu+1)+cu*self.beta(2+au+0.5,bu+1)+du*self.beta(2+au+1.0,bu+1)
    normd = self.beta(2+ad,bd+1)+cd*self.beta(2+ad+0.5,bd+1)+dd*self.beta(2+ad+1.0,bd+1)
    uv = lambda y: Nu*y**au*(1-y)**bu*(1+cu*y**0.5+du*y)/normu
    dv = lambda y: Nd*y**ad*(1-y)**bd*(1+cd*y**0.5+dd*y)/normd
    qv = lambda y: uv(y) + dv(y)
    integrand = lambda y: func(y) * qv(y)
    sumrule = quad(integrand,0,1.0)[0]
    print 'The sum rule is %5.10f (direct PDF form)'%sumrule

    #--test that sum rule is satisfied with inverse Mellin transform of pdf

    _qv = lambda y: (self.pdf.get_xF(y,Q2,'uv')+self.pdf.get_xF(y,Q2,'dv'))/y
    integrand = lambda y: func(y) * _qv(y)
    sumrule = quad(integrand,0,1.0)[0]
    print 'The sum rule is %5.10f (Mellin transform of PDF)'%sumrule

    #--test that function matches inverse mellin

    print 'xF = %5.6f (direct form)'%(x*func(x))
    print 'xF = %5.6f (Mellin transform)'%( x*self.mellin.invert(x,self.off[nucleon]))

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

    offshell.test(0.1,1.27**2,'proton')

    offshell.test(0.1,10,'proton')





