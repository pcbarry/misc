#!/usr/bin/env python
from mpmath import fp
import numpy as np
from scipy.integrate import quad
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
    offshell(x) = norm*x^a*(1-x)^b*(1+c*x^0.5+d*x)
    """
    params={}
    params['proton']  =np.array([1.0,-0.5,10.0,1.0,1.0])
    params['neutron'] =np.array([1.0,0.5,3.0,1.0,1.0])
    self.params=params

  def set_sumrules(self):
    vmom={}
    for flav in self.params:
        vmom[flav]={}
        a,b,c=self.params[flav][1],self.params[flav][2],self.params[flav][3]
        vmom[flav]['a+1']=self.pdf.get_moments('uv1',N=a+1,b_shift=b)+self.pdf.get_moments('dv1',N=a+1,b_shift=b)
        vmom[flav]['a+1.5']=self.pdf.get_moments('uv1',N=a+1.5,b_shift=b)+self.pdf.get_moments('dv1',N=a+1.5,b_shift=b)
        vmom[flav]['a+2']=self.pdf.get_moments('uv1',N=a+2,b_shift=b)+self.pdf.get_moments('dv1',N=a+2,b_shift=b)

        d=(-vmom[flav]['a+1']-c*vmom[flav]['a+1.5'])/(vmom[flav]['a+2'])
        self.params[flav][4]=d

  def get_moments(self,flav,M=None):
    """
    if M==None: then parametrization is to be use to compute moments along mellin contour
    else the Mth moment is returned
    """
    if np.array(M)[0]==None: M=self.dmellin.M
    norm,a,b,c,d=self.params[flav]
    return norm*(self.aux.beta(M+a,b+1)+c*self.aux.beta(M+a+0.5,b+1)+d*self.aux.beta(M+a+1,b+1))

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

  def test(self,x,Q2,nucleon):
    norm,a,b,c,d=self.params[nucleon]
    print d
    func = lambda y: norm*(y**a)*((1-y)**b)*(1+c*y**0.5+d*y)
    self.pdf.setup()
    self.pdf.evolve(Q2)

    print 'At Q2 = %5.4f:'%Q2

    #--test that sum rule is satisfied with parameterized pdf form

    Nu,au,bu,cu,du=self.pdf.params['uv1']
    Nd,ad,bd,cd,dd=self.pdf.params['dv1']
    uv = lambda y: Nu*y**au*(1-y)**bu
    dv = lambda y: Nd*y**ad*(1-y)**bd
    qv = lambda y: uv(y) + dv(y)
    integrand = lambda y: func(y) * qv(y)
    sumrule = quad(integrand,0,1.0)[0]
    print 'The sum rule is %5.10f (direct PDF form)'%sumrule

    #--test that sum rule is satisfied with inverse Mellin transform of pdf

    qv = lambda y: (self.pdf.get_xF(y,Q2,'uv')+self.pdf.get_xF(y,Q2,'dv'))/y
    #qv = lambda y: (self.pdf.get_xF0(y,'um')+self.pdf.get_xF0(y,'dm'))/y
    integrand = lambda y: func(y) * qv(y)
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





