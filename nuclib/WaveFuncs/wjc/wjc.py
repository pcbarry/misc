#!/usr/bin/env python
import sys,os
import numpy as np
from scipy.interpolate import griddata

#--deuteron wavefunctions from WJC (Gross et al.) NN potential model
#--Ref: Gross & Stadler, Phys. Rev. C78, 014005 (2008)

#--from data: q in MeV, wavefunctions in GeV^-3/2
#--ouput:     q in 1/fm, wavefunctions in fm^3/2

#--Wave function normalization: integral(dq q^2 (u^2 + w^2 +vt^2 + vs^2) + V' term) = 1
#--Wave function is normalized to ~105.479% for wjc-1
#--Wave function is normalized to ~102.31% for wjc-2
#--For structure functions renormalize to 1

def wjc(k,kind=1):
    hcM = 197.327  #--conversion factor MeV*fm
    hcG = 0.197327 #--conversion factor GeV*fm

    #--load data
    if kind==1: data = np.load('WaveFuncs/wjc/wjc-1.npy')
    if kind==2: data = np.load('WaveFuncs/wjc/wjc-2.npy')

    k_data  = data[0]/hcM       #--convert from MeV to 1/fm
    u_data  = data[1]*hcG**1.5  #--convert from GeV^-3/2 to fm^3/2
    w_data  = data[2]*hcG**1.5

    ##--test normalization 
    #x,w = np.polynomial.legendre.leggauss(99)
    #kmin = np.min(k_data)
    ##--important contributions to integral go up to ~80
    #kmax = 100

    #_k  = 0.5*(x*(kmax - kmin) + kmin + kmax) 
    #jac = 0.5*(kmax - kmin)

    #_u = lambda k: griddata(k_data,u_data,k,method='cubic',fill_value=0)
    #_w = lambda k: griddata(k_data,w_data,k,method='cubic',fill_value=0)
    #_vt = lambda k: griddata(k_data,vt_data,k,method='cubic',fill_value=0)
    #_vs = lambda k: griddata(k_data,vs_data,k,method='cubic',fill_value=0)
  
    #_integrand = lambda k: k**2 * ((_u(k))**2 + (_w(k))**2 + (_vt(k))**2 + (_vs(k))**2) 

    #integrand = np.zeros(len(x))
    #for i in range(len(x)):
    #    integrand[i] = _integrand(_k[i]) 

    #N = np.sum(w*integrand*jac)

    #print N
 
    u = griddata(k_data,u_data,k,method='cubic',fill_value=0) 
    w = griddata(k_data,w_data,k,method='cubic',fill_value=0) 

    return u,w

















