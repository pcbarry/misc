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
    offshell(x) = norm*(x-x0)(x-x1)(1+x0-x)
    No sum rule!
    """
    params={}
    params['proton']  =np.array([1,0.5,3])
    params['neutron'] =np.array([1,0.5,3])
    self.params=params

    vmom={}
    for flav in self.params:
        vmom[flav]={}
        x0=self.params[flav][1]
        for i in range(4):
            vmom[flav][i+1]=self.pdf.get_moments('uv1',i+1)+self.pdf.get_moments('dv1',i+1)

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

  def get_xF(self,x,nucleon):
    return x*self.mellin.invert(x,self.off[nucleon])

if __name__=='__main__':

    conf['order']='NLO'
    conf['Q20'] = 1.0
    conf['aux']=AUX()
    conf['dglap mode']='truncated'
    conf['mellin']=MELLIN(npts=8)
    conf['dmellin']=DMELLIN(nptsN=8,nptsM=8)
    conf['alphaS']=ALPHAS()
    conf['pdf']=pdf2.PDF()

    offshell=OFFSHELL()

    offshell.test(0.1,1.27**2,'proton')






