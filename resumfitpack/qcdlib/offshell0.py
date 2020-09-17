#!/usr/bin/env python
from mpmath import fp
import numpy as np
from qcdlib.alphaS import ALPHAS
from qcdlib.aux import AUX
from qcdlib.mellin import MELLIN
from tools.config import conf

  
class OFFSHELL: 

  def __init__(self):
    self.aux=conf['aux']
    self.mellin=conf['mellin']
    self.dmellin=conf['dmellin']
    self.set_params()
    self.setup()

  def set_params(self):
    params={}
    params['proton']  =np.array([1,0.5,3,0,0])
    params['neutron'] =np.array([1,0.5,3,0,0])
    self.params=params

  def get_moments(self,flav,M=None):
    """
    if M==None: then parametrization is to be use to compute moments along mellin contour
    else the Mth moment is returned
    """
    if np.array(M)[0]==None: N=self.dmellin.M
    norm,a,b,c,d=self.params[flav]
    return norm*(self.aux.beta(N+a,b+1)+c*self.aux.beta(N+a+0.5,b+1)+d*self.aux.beta(N+a+1,b+1))

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
    conf['mellin']=MELLIN(npts=8)
    conf['alphaS']=ALPHAS()

    t4=T4()

    print OFFSHELL.get_xF(0.5,10.,'proton')






