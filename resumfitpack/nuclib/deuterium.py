#!/usr/bin/env python
import sys
import numpy as np
from tools.tools import checkdir,load,save
from scipy.interpolate import RectBivariateSpline,InterpolatedUnivariateSpline
from tools.bar import BAR
from scipy.integrate import quad,fixed_quad
from scipy.optimize import root
from multiprocessing import Process,Queue,Pool,Pipe
from scipy.interpolate import griddata 
try:import vegas
except: pass
try: from WaveFuncs.paris.paris import paris 
except: pass
import time

class DEUTERON: 
    """
    Notes:  
      - y range is [0,mD/mN]
      - we need to normalize f0 such that int_0^2 dy f22(y,gam=1) = 1    
    """
  
    def __init__(self,path2nuctab='.',precalc=False):
  
        self.path2nuctab=path2nuctab.rstrip('/')
        self.D = {}
        self.setup_internal()
        if precalc==False: 
            self.load_tables()
  
    # intermediate steps
  
    def setup_internal(self):
        self.hc = 197.327                  # conversion factor (MeV.fm)
        self.mN = 938.91897                # nucleon mass
        self.epsD = -2.224575              # deuteron binding energy
        self.mD = 2*self.mN + self.epsD    # deuteron mass
        self.mN2=self.mN**2
        self.mD2=self.mD**2
        self.pTmin=0
        self.pTmax=20000
        self.mNfm=self.mN/self.hc
        self.epsDfm=self.epsD/self.hc
        self.mN2_=0.93891897**2            # for interpolation
        self.ymin=0
        self.ymax=self.mD/self.mN
        self.norm=1.0 
  
    def get_kin(self,y,gam,pT):
      
        yD = y*self.mN/self.mD   # Convert to y -> yD
  
        # Use relativistic kinematics 
        if gam>1:
            pz = (np.sqrt((1.-yD)**2*self.mD2+(gam**2-1.)*(pT**2+self.mN2)) - self.mD*(1.-yD)*gam )/(gam**2-1.)
        else:
            y0 = yD
            pz = (pT**2 + self.mN2 - (1.-y0)**2*self.mD2) / (2.*(1-y0)*self.mD)
  
        pv  = np.sqrt(pT**2 + pz**2)   # Total nucleon 3-momentum
        pv2 = pT**2 + pz**2
        Ep  = np.sqrt(self.mN2 + pv2)  # On-shell spectator nucleon energy
        p0  = self.mD - Ep             # Off-shell interacting nucleon energy
        p2  = p0**2 - pv2              # Nucleon virtuality
        pT2 = pT**2
        eps = self.mD - Ep - self.mN   # Separation energy
  
        # Bjorken limit light-cone variable
        # Solution of y = y0/2 (1 + gamma - (gamma-1)/y0^2 (p2+pT2)/M^2)
        y0 = (yD+np.sqrt(yD**2 +(gam**2-1)*(p2+pT**2)/self.mD2))/(1.+gam) 
  
        cos = np.sqrt(pz**2/pv2)   # Cosine squared (pz_hat in formula)
  
        return pv,p2,pz,pT2,cos,y0,Ep,eps
  
    def get_upol_Dij(self,y,gam,pT,fXX,kind,f0):
        p,p2,pz,pT2,cos,y0,Ep,eps = self.get_kin(y,gam,pT)
  
        if   'fLL' in fXX:  C_ij = 1.
        elif 'fL2' in fXX:  
            if gam==1:  C_ij = 0.0
            else:       C_ij = (gam**2-1.)*pT2/(y*self.mN)**2
        elif 'f22' in fXX:
            if gam==1:  C_ij = 1.
            else:       C_ij = 1./gam**2*(1.+(gam**2-1.)/(2.*y**2*self.mN**2)*(2.*p2+3.*pT2))
        elif 'f33' in fXX:  C_ij = 1.
        else:  raise ValueError('fXX not available at get_upol_Dij')
  
        if 'f33' in fXX: flux = 1. + pz/(gam*self.mN) 
        else:            flux = 1. + pz*gam/self.mN
  
        Dij=flux*C_ij*f0
  
        if kind=='onshell':  Dij*=1
        elif kind=='offshell': Dij*=(p2-self.mN2)/self.mN2
  
        return Dij
  
    def get_pol_Dij(self,y,gam,pT,fXX,kind,f1,f2):
        p,p2,pz,pT2,cos,y0,Ep,eps = self.get_kin(y,gam,pT)
  
        if fXX=='g11':
            Dij =   f1+(3.-gam**2)/(6.*gam**2)*(3.*cos**2-1.)\
                    *f2+(pz/gam/self.mN)*(f1+2.*f2/3.)+(p/self.mN)**2\
                    *((3.-gam**2)*cos**2-1.-gam**2)/(12.*gam**2)\
                    *(3.*f1-f2)
        elif fXX=='g12':
            Dij =   (gam**2-1.)*((1.-3.*cos**2)/(2.*gam**2)*f2\
                    +(pz/gam/self.mN)*(f1+(3.*cos**2/2.-5./6.)*f2)\
                    -(p/self.mN)**2*((1.+cos**2*(4.*gam**2-3.))/(4.*gam**2)\
                    *f1+(5.+18.*cos**4*gam**2-5.*cos**2\
                    *(3.+2.*gam**2))/(12.*gam**2)*f2))
        elif fXX=='g21':
            Dij =   (1.-3.*cos**2)/(2.*gam**2)*f2-(pz/gam/self.mN)\
                    *(f1+2.*f2/3.)-(p/self.mN)**2*(3.*cos**2-1.)\
                    /(12.*gam**2)*(3.*f1-f2)
        elif fXX=='g22':
            Dij =   f1+(2.*gam**2-3.)/(6.*gam**2)\
                    *(3.*cos**2-1.)*f2+(pz/gam/self.mN)*((1.-gam**2)\
                    *f1+(gam**2/3.0-5./6.+cos**2*(3./2.-gam**2))\
                    *f2)+(p/self.mN)**2*((cos**2*(3.-6.*gam**2+4.*gam**4)-1.\
                    -2.*gam**2)/(4.*gam**2)*f1+(5.-2.*gam**2\
                    *(1.+3.*cos**2)+4.*cos**2*gam**4)/(12.*gam**2)\
                    *(3.*cos**2-1.)*f2)
        else:
            raise ValueError('fXX not available at get_pol_Dij')
  
        if kind=='onshell':  Dij*=1
        elif kind=='offshell': Dij*=(p2-self.mN2)/self.mN2
  
        return Dij
  
    def get_y0(self,y,gam,pT):
        p,p2,pz,pT2,cos,y0,Ep,eps = self.get_kin(y,gam,pT)
        return y0
  
    def integrand(self,X,gam,fXX,kind):
  
        y=self.ymin+X[1]*(self.ymax-self.ymin)
  
        if gam==1:
            pTmax=self.pTmax
        else:
            sol = root(lambda _pT: 1-self.get_y0(y,gam,_pT), [self.pTmax],  method='hybr')
            pTmax=sol.x[0]
  
        pT=self.pTmin+X[0]*(pTmax-self.pTmin)
  
        # get kinematics
        p,p2,pz,pT2,cos,y0,Ep,eps = self.get_kin(y,gam,pT)
  
        # Deuteron wave function (Momenta in MeV)
        # Nonrelativistic wave functions 
        u,w =paris(1,p/self.hc)
  
        # Spectral Function Coefficients
        # Output wavefunctions in fm^3/2 => MeV^-3/2
        u/= self.hc**1.5
        w/= self.hc**1.5
        f0 = (u**2 + w**2)/self.norm
        f1 = (u**2 - w**2/2.)
        f2 = (3./2.)*(w**2 - np.sqrt(2.)*w*u)
  
        if   fXX.startswith('f'):   Dij=self.get_upol_Dij(y,gam,pT,fXX,kind,f0)
        elif fXX.startswith('g'):   Dij=self.get_pol_Dij(y,gam,pT,fXX,kind,f1,f2)   
        else:  raise ValueError('fXX not available')
  
        jac = pT*Ep*self.mN/gam/self.mD/(1.-y0)
        return 0.5*jac*Dij*(pTmax-self.pTmin)
  
    def y2x(self,y):
        return (y-self.ymin)/(self.ymax-self.ymin)
  
    def set_normalization(self):
        print 'computing normalization...'
  
        gam=1
        kind='onshell'
        fXX='f22'
  
        df10 = np.vectorize(lambda x1,x0:self.integrand([x0,x1],gam,fXX,kind))
        df1  = np.vectorize(lambda x1: fixed_quad(lambda x0: df10(x1,x0),0,1,n=1000)[0])
        self.norm = fixed_quad(lambda x1: df1(x1),0,1,n=200)[0]*(self.ymax-self.ymin)
        print 'norm=',self.norm
  
    def _get_fXX(self,y,gam,fXX,kind):
        x1=self.y2x(y)
        f=0
        df0 = np.vectorize(lambda x0:self.integrand([x0,x1],gam,fXX,kind))
        return fixed_quad(lambda x0: df0(x0),0,1,n=1000)[0]
  
    # subroutines for table generation and interpolation 
  
    def gen_table(self,G,Y,fXX,kind):
  
        bar=BAR('computing %s-%s '%(fXX,kind),G.size*Y.size)
        table = np.zeros((G.size,Y.size))
        for iG in range(G.size):
          for iY in range(Y.size): 
            table[iG,iY]=self._get_fXX(Y[iY],G[iG],fXX,kind)
            bar.next()
        bar.finish()
        save(table,'%s/%s-%s.dat'%(self.path2nuctab,fXX,kind))
  
    def gen_tables(self):
  
        checkdir(self.path2nuctab)
        G = np.linspace(1.,4.,31)
        Y = np.linspace(0.01,1.99,199)
        save(Y,'%s/Y.dat'%(self.path2nuctab))
        save(G,'%s/G.dat'%(self.path2nuctab))
  
        t1=time.time()
        #P = [Process(target=self.gen_table, args=(G,Y,fXX,'onshell')) for fXX in ['f22','fLL','fL2','f33']]
        #P.extend([Process(target=self.gen_table, args=(G,Y,fXX,'offshell')) for fXX in ['f22','fLL','fL2','f33']])
        P = [Process(target=self.gen_table, args=(G,Y,fXX,'onshell')) for fXX in ['g11','g12','g21','g22']]
        P.extend([Process(target=self.gen_table, args=(G,Y,fXX,'offshell')) for fXX in ['g11','g12','g21','g22']])
        for p in P: p.start() 
        for p in P: p.join() 
        t2=time.time()
  
    def load_tables(self):
  
        Y = load('%s/Y.dat'%(self.path2nuctab))
        G = load('%s/G.dat'%(self.path2nuctab))
        self.Y,self.G=[],[]
        for iG in range(G.size):
            for iY in range(Y.size): 
                self.Y.append(Y[iY])
                self.G.append(G[iG])
        self.Y=np.array(self.Y)
        self.G=np.array(self.G)

        self.tabs={}
  
        fXXs=['f22','fLL','fL2','f33','g11','g12','g21','g22']
        for fXX in fXXs:
            for kind in ['onshell','offshell']:
                tab=load('%s/%s-%s.dat'%(self.path2nuctab,fXX,kind))
                self.D['%s-%s'%(fXX,kind)]=RectBivariateSpline(G,Y,tab,kx=1, ky=1)
                self.tabs[(fXX,kind)]=[]
                for iG in range(G.size):
                    for iY in range(Y.size): 
                        self.tabs[(fXX,kind)].append(tab[iG,iY])
                self.tabs[(fXX,kind)]=np.array(self.tabs[(fXX,kind)])

    def get_fXX(self,fXX,kind,x,Q2,y):
        g=(1+4*x**2*self.mN2_/Q2)**0.5
        return self.D['%s-%s'%(fXX,kind)](g,y)[0,0]

    def get_fXX2(self,fXX,kind,x,Q2,y):
        g=np.array((1+4*x**2*self.mN2_/Q2)**0.5)
        tab=self.tabs[(fXX,kind)]
        return griddata((self.G,self.Y),tab,(g ,y),fill_value=0, method='cubic')  


  
if __name__== "__main__":

    deuteron=DEUTERON(path2nuctab='./grids/deuteron',precalc=True)
    #deuteron.set_normalization()
    #deuteron.gen_tables()
    
    deuteron.load_tables()
    print(deuteron.get_fXX('f22','onshell',0.5,10.,0.6))


