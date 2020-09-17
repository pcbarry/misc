#!/usr/bin/env python
import sys
import numpy as np
from tools.tools import checkdir,load,save
from scipy.interpolate import RectBivariateSpline,InterpolatedUnivariateSpline, griddata
from tools.bar import BAR
from scipy.integrate import quad,fixed_quad
from multiprocessing import Process,Queue,Pool,Pipe
try:import vegas
except: pass
import time

class HELIUM:
  
    def __init__(self,path2nuclib='.',path2nuctab='.',group='kpsv',precalc=False):
  
        self.group=group
        self.path2nuclib=path2nuclib.rstrip('/')
        self.path2nuctab=path2nuctab.rstrip('/')
        self.D = {}
        self.setup_internal()
  
        if precalc:                        # Load 3He spectral function data
            self.load_table2D('f0pc')
            self.load_table2D('f1pc')
            self.load_table2D('f2pc')
            self.load_table2D('f0nc')
            self.load_table2D('f1nc')
            self.load_table2D('f2nc')
            self.load_table1D('f0pd')
            self.load_table1D('f1pd')
            self.load_table1D('f2pd')
        else: 
            self.load_tables()
  
    # routines to load spectral functions
  
    def load_table2D(self,name):
        '''Load spectral function data (cont).'''
        '''Must be in units of MeV!'''   
    
        F = open('%s/WaveFuncs/%s/%s.dat'%(self.path2nuclib,self.group,name),'r')
        L = F.readlines()
        F.close()
    
        L = [l.strip() for l in L]
        L = [[float(x) for x in l.split()] for l in L]
        L = np.transpose(L)
  
        EPS = np.array(sorted(set(L[0])))  # what is this?
        MOM = np.array(sorted(set(L[1])))  # what is this?
  
        SF_ = L[2]
        SF = np.zeros((EPS.size,MOM.size))
  
        cnt = 0
        for iEPS in range(EPS.size):
            for iMOM in range(MOM.size):
                SF[iEPS,iMOM] = SF_[cnt]
                cnt += 1
  
        self.D[name] = RectBivariateSpline(EPS,MOM,SF)
  
        self.emn = EPS[0]
        self.emx = EPS[-1]
        self.pmn = MOM[0]
        self.pmx = MOM[-1]
  
    def load_table1D(self,name):
        '''Load spectral function data (pole).'''
        '''Must be in units of MeV!'''   
 
        F = open('%s/WaveFuncs/%s/%s.dat'%(self.path2nuclib,self.group,name),'r')
        L = F.readlines()
        F.close()
    
        L = [l.strip() for l in L]
        L = [[float(x) for x in l.split()] for l in L]
        L = np.transpose(L)
  
        MOM = L[1]
        SF = L[2]

        #--files f0pd, f1pd, f2pd in WaveFuncs/kpsv need to be converted from 1/fm to MeV 
        if self.group=='kpsv':
            convert = ['f0pd','f1pd','f2pd']
            if name in convert:
                SF*=self.hc
 
        self.D[name] = InterpolatedUnivariateSpline(MOM,SF)
        self.pmn = MOM[0]
        self.pmx = MOM[-1]
  
    # routines to compute fXX
  
    def setup_internal(self):
        # some explanation  for each constance seems ideal !!!!
        self.hc = 197.327053   # 
        self.mN = (938.27231 + 939.56563)/2.0/self.hc
        self.mN2 = self.mN**2
        self.eD = -2.2246/self.hc
        self.mD = 2.0*self.mN + self.eD
        self.eHe3 = -7.72/self.hc
        self.m0 = 2.0*self.mN
        self.e0p = self.eHe3 - self.eD
        self.e0c = self.eHe3
  
        self.ymin=0
        self.ymax=3.0
  
        self.norm_p = 1
        self.norm_n = 1
        self.a=-4.9460565      # cont e integration mapping 
  
        # for interpolation
        self.mN2_=0.93891897**2   # For the calculation mN is used. For interpolation mN_ is used. is this right? 
    
    def get_kin(self,y,gam,p,eps):
        cos = (self.mN * (y - 1.0) - eps)/(gam*p)
        pz   = p*cos
        p2   = (self.mN+eps)**2 - p**2
        pT2  = p**2-pz**2
        return p,p2,pz,pT2,cos
  
    def get_upol_Dij(self,y,gam,p,eps,fXX,kind,f0):
        p,p2,pz,pT2,cos = self.get_kin(y,gam,p,eps)
  
        if   'fLL' in fXX:  C_ij = 1.
        elif 'fL2' in fXX:  
            if gam==1:  C_ij = 0.0
            else:       C_ij = (gam**2-1.)*pT2/(y*self.mN)**2
        elif 'f22' in fXX:
            if gam==1:  C_ij = 1.
            else:       C_ij = 1./gam**2*(1.+(gam**2-1.)/(2.*y**2*self.mN**2)*(2.*p2+3.*pT2))
        elif 'f3'  in fXX:  C_ij = 1.
        else:  raise ValueError('fXX not available at get_upol_Dij')
  
        if 'f3' in fXX:  flux = 1. + pz/(gam*self.mN)   
        else:            flux = 1. + pz*gam/self.mN
  
        Dij=flux*C_ij*f0
  
        if kind=='onshell':  Dij*=1
        elif kind=='offshell': Dij*=(p2-self.mN2)/self.mN2
  
        return Dij
  
    def get_pol_Dij(self,y,gam,p,eps,fXX,kind,f1,f2):
        p,p2,pz,pT2,cos = self.get_kin(y,gam,p,eps)
        mN = self.mN 
 
        if 'g11' in fXX:
            Dij =   f1+(3.-gam**2)/(6.*gam**2)*(3.*cos**2-1.)\
                    *f2+(pz/gam/mN)*(f1+2.*f2/3.)+(p/mN)**2\
                    *((3.-gam**2)*cos**2-1.-gam**2)/(12.*gam**2)\
                    *(3.*f1-f2)
        elif 'g12' in fXX:
            Dij =   (gam**2-1.)*((1.-3.*cos**2)/(2.*gam**2)*f2\
                    +(pz/gam/mN)*(f1+(3.*cos**2/2.-5./6.)*f2)\
                    -(p/mN)**2*((1.+cos**2*(4.*gam**2-3.))/(4.*gam**2)\
                    *f1+(5.+18.*cos**4*gam**2-5.*cos**2\
                    *(3.+2.*gam**2))/(12.*gam**2)*f2))
        elif 'g21' in fXX:
            Dij =   (1.-3.*cos**2)/(2.*gam**2)*f2-(pz/gam/mN)\
                    *(f1+2.*f2/3.)-(p/mN)**2*(3.*cos**2-1.)\
                    /(12.*gam**2)*(3.*f1-f2)
        elif 'g22' in fXX:
            Dij =   f1+(2.*gam**2-3.)/(6.*gam**2)\
                    *(3.*cos**2-1.)*f2+(pz/gam/mN)*((1.-gam**2)\
                    *f1+(gam**2/3.0-5./6.+cos**2*(3./2.-gam**2))\
                    *f2)+(p/mN)**2*((cos**2*(3.-6.*gam**2+4.*gam**4)-1.\
                    -2.*gam**2)/(4.*gam**2)*f1+(5.-2.*gam**2\
                    *(1.+3.*cos**2)+4.*cos**2*gam**4)/(12.*gam**2)\
                    *(3.*cos**2-1.)*f2)
        else:  raise ValueError('fXX not available at get_pol_Dij')
  
        if kind=='onshell':  Dij*=1
        elif kind=='offshell': Dij*=(p2-self.mN2)/self.mN2
  
        return Dij
  
    def get_plims(self,a,gam):
        z1 = 8.0*gam**2-4.0*a-np.sqrt((8.0*gam**2-4.0*a)**2-(4.0*a)**2)
        z2 = 8.0*gam**2-4.0*a+np.sqrt((8.0*gam**2-4.0*a)**2-(4.0*a)**2)
        pmin = np.amax([self.pmn,self.mN*np.sqrt(z1)])
        pmax = np.amin([self.pmx,self.mN*np.sqrt(z2)])
        return pmin,pmax
  
    def _integrand_pole(self,y,gam,fXX,kind,p):
        '''
        recall that this is only for fXXp
        '''
        eps = self.eHe3 - self.eD - p**2/(2.0*self.mD)
  
        if   fXX.startswith('f'):
  
            f0 = 4.0*np.pi/self.hc*self.D['f0pd'](p)/self.norm_p
            Dij=self.get_upol_Dij(y,gam,p,eps,fXX,kind,f0) 
  
        elif fXX.startswith('g'):
  
            f1 = 4.0*np.pi/self.hc*self.D['f1pd'](p)
            f2 = 4.0*np.pi/self.hc*self.D['f2pd'](p)
            Dij=self.get_pol_Dij(y,gam,p,eps,fXX,kind,f1,f2)   
  
        else:
            raise ValueError('fXX not available at integrand_pole')
  
        jac  = self.mN*p/gam
  
        return 0.5*jac*Dij
  
    def _integrand_cont(self,y,gam,fXX,kind,p,e):
  
        eps = self.eHe3 - e/self.hc - p**2/(2*self.m0)
  
        if   fXX.startswith('f'):
  
            if fXX.endswith('p'): f0 = 4.0*np.pi/self.hc*self.D['f0pc'](e,p)[0,0]/self.norm_p
            if fXX.endswith('n'): f0 = 4.0*np.pi/self.hc*self.D['f0nc'](e,p)[0,0]/self.norm_n
            Dij=self.get_upol_Dij(y,gam,p,eps,fXX,kind,f0) 
  
        elif fXX.startswith('g'):
  
            if fXX.endswith('p'):
                f1 = 4.0*np.pi/self.hc*self.D['f1pc'](e,p)[0,0]
                f2 = 4.0*np.pi/self.hc*self.D['f2pc'](e,p)[0,0]
            if fXX.endswith('n'):
                f1 = 4.0*np.pi/self.hc*self.D['f1nc'](e,p)[0,0]
                f2 = 4.0*np.pi/self.hc*self.D['f2nc'](e,p)[0,0]
            Dij=self.get_pol_Dij(y,gam,p,eps,fXX,kind,f1,f2)   
  
        else:
            raise ValueError('fXX not available at integrand_pole')
  
        jac  = self.mN*p/gam
        return 0.5*jac*Dij
  
    def integrand_pole(self,X,gam,fXX,kind):
        y=self.ymin+X[1]*(self.ymax-self.ymin)
        a = y - 1.0 - self.e0p/self.mN
        if a<=gam**2:
            pmin,pmax=self.get_plims(a,gam)
            p=pmin + X[0]*(pmax-pmin)
            jac=(pmax-pmin)
            return jac* self._integrand_pole(y,gam,fXX,kind,p)
        else:
            return 0
  
    def integrand_cont(self,X,gam,fXX,kind):
        y=self.ymin+X[2]*(self.ymax-self.ymin)
        elim = self.hc*(self.mN*(gam**2+1.0-y)+self.e0c)
        if elim>self.emn:
            emaxc = np.amin([elim,self.emx])  
            eminc = self.emn
            e=eminc + X[1]*(emaxc-eminc)
            emev = -self.e0c * self.hc + e
            a = y - 1.0 + (emev/self.hc)/self.mN
            pmin,pmax=self.get_plims(a,gam)
            p=pmin + X[0]*(pmax-pmin)
            jac= (pmax-pmin) * (emaxc-eminc)  
            return jac * self._integrand_cont(y,gam,fXX,kind,p,e)
        else:
            return 0
  
    def y2x(self,y):
        return (y-self.ymin)/(self.ymax-self.ymin)
  
    def set_normalizations(self):
        print 'computing normalizations...'
  
        gam=1
        kind='onshell'  # normalizations are set by the onshell contribution
  
        for fXX in ['f22p','f22n']:
  
            if fXX=='f22p': 
                dfp10 = np.vectorize(lambda x1,x0:self.integrand_pole([x0,x1],gam,fXX,kind))
                dfp1  = np.vectorize(lambda x1: fixed_quad(lambda x0: dfp10(x1,x0),0,1,n=100)[0])
                dfp   = fixed_quad(lambda x1: dfp1(x1),0,1,n=100)[0]*(self.ymax-self.ymin)
  
            dfc210 = np.vectorize(lambda x2,x1,x0:self.integrand_cont([x0,x1,x2],gam,fXX,kind))
            dfc21  = np.vectorize(lambda x2,x1: fixed_quad(lambda x0: dfc210(x2,x1,x0),0,1,n=100)[0])
            dfc2   = np.vectorize(lambda x2: fixed_quad(lambda x1: dfc21(x2,x1),0,1,n=200)[0])
            dfc    = fixed_quad(lambda x2: dfc2(x2),0,1,n=50)[0]*(self.ymax-self.ymin)

            if fXX=='f22p': 
                self.norm_p=dfp+dfc
                print 'norm_p=',self.norm_p
            if fXX=='f22n': 
                self.norm_n=dfc
                print 'norm_n=',self.norm_n
  
    def _get_fXX(self,y,gam,fXX,kind,part='all'):
        x=self.y2x(y)
        f=0
  
        if part=='all' or part=='pole':
            # compute pole contribution 
            if fXX.endswith('p'):
                dfp0 = np.vectorize(lambda x0:self.integrand_pole([x0,x],gam,fXX,kind))
                f+= fixed_quad(lambda x0: dfp0(x0),0,1,n=100)[0]
  
        if part=='all' or part=='cont':        
            # compute cont contribution 
            dfc10 = np.vectorize(lambda x1,x0:self.integrand_cont([x0,x1,x],gam,fXX,kind))
            dfc1  = np.vectorize(lambda x1: fixed_quad(lambda x0: dfc10(x1,x0),0,1,n=100)[0])
            f+= fixed_quad(lambda x1: dfc1(x1),0,1,n=200)[0]
  
        return f
  
    # routines for fXX table generation and interpolation
  
    def _gen_table(self,G,Y,fXX,kind,part):
        bar=BAR('computing %s-%s '%(fXX,self.group),G.size*Y.size)
        table = np.zeros((G.size,Y.size))
        for iG in range(G.size):
            for iY in range(Y.size): 
                table[iG,iY]=self._get_fXX(Y[iY],G[iG],fXX,kind,part)
                bar.next()
        bar.finish()
        if part=='all':
            save(table,'%s/%s-%s-%s.dat'%(self.path2nuctab,fXX,kind,self.group))
        else:
            save(table,'%s/%s-%s-%s-%s.dat'%(self.path2nuctab,fXX,kind,part,self.group))
  
    def gen_tables(self): 
        print 'computing tables...'
        checkdir(self.path2nuctab)
        G = np.linspace(1.,4.,31)
        Y = np.linspace(0.03,self.ymax,300)
        save(Y,'%s/Y.dat'%(self.path2nuctab))
        save(G,'%s/G.dat'%(self.path2nuctab))
  
        #fXXs=['f22p','fLLp','fL2p','f22n','fLLn','fL2n']
        #P = [Process(target=self._gen_table, args=(G,Y,fXX,'onshell','all',)) for fXX in fXXs]
        #P.extend([Process(target=self._gen_table, args=(G,Y,fXX,'offshell','all',)) for fXX in fXXs])
        #for p in P: p.start() 
        #for p in P: p.join() 
  
  
        #fXXs=['f22p','fLLp','fL2p']
        #P = [Process(target=self._gen_table, args=(G,Y,fXX,'onshell','pole')) for fXX in fXXs]
        #P.extend([Process(target=self._gen_table, args=(G,Y,fXX,'offshell','pole',)) for fXX in fXXs])
        #for p in P: p.start() 
        #for p in P: p.join() 
  
        #fXXs=['f22p','fLLp','fL2p','f22n','fLLn','fL2n']
        #P = [Process(target=self._gen_table, args=(G,Y,fXX,'onshell','cont')) for fXX in fXXs]
        #P.extend([Process(target=self._gen_table, args=(G,Y,fXX,'offshell','cont',)) for fXX in fXXs])
        #for p in P: p.start() 
        #for p in P: p.join()
    
        #--polarized structure functions 
        fXXs=['g11p','g12p','g21p','g22p','g11n','g12n','g21n','g22n']
        P = [Process(target=self._gen_table, args=(G,Y,fXX,'onshell','all',)) for fXX in fXXs]
        P.extend([Process(target=self._gen_table, args=(G,Y,fXX,'offshell','all',)) for fXX in fXXs])
        for p in P: p.start() 
        for p in P: p.join() 
  
        #fXXs=['g11p','g12p','g21p','g22p']
        #P = [Process(target=self._gen_table, args=(G,Y,fXX,'onshell','pole')) for fXX in fXXs]
        #P.extend([Process(target=self._gen_table, args=(G,Y,fXX,'offshell','pole',)) for fXX in fXXs])
        #for p in P: p.start() 
        #for p in P: p.join() 
  
        #fXXs=['g11p','g12p','g21p','g22p','g11n','g12n','g21n','g22n']
        #P = [Process(target=self._gen_table, args=(G,Y,fXX,'onshell','cont')) for fXX in fXXs]
        #P.extend([Process(target=self._gen_table, args=(G,Y,fXX,'offshell','cont',)) for fXX in fXXs])
        #for p in P: p.start() 
        #for p in P: p.join() 
  
    def load_tables(self):
 
        Y = load('%s/Y.dat'%(self.path2nuctab))
        G = load('%s/G.dat'%(self.path2nuctab))
 
        self.Y,self.G=[],[]
        for iG in range(G.size):
            for iY in range(Y.size):
                self.Y.append(Y[iY])
                self.G.append(G[iG])
        self.Y = np.array(self.Y)
        self.G = np.array(self.G)

        self.tabs = {}
 
        fXXs=['f22p','fLLp','fL2p','f22n','fLLn','fL2n']
        fXXs.extend(['g11p','g12p','g21p','g22p','g11n','g12n','g21n','g22n'])
        for fXX in fXXs:
            for kind in ['onshell','offshell']:
                tab=load('%s/%s-%s-%s.dat'%(self.path2nuctab,fXX,kind,self.group))
                self.D['%s-%s'%(fXX,kind)]=RectBivariateSpline(G,Y,tab)
                self.tabs[(fXX,kind)] = []
                for iG in range(G.size):
                    for iY in range(Y.size):
                        self.tabs[(fXX,kind)].append(tab[iG,iY])
                self.tabs[(fXX,kind)] = np.array(self.tabs[(fXX,kind)]) 
 
        fXXs=['f22p','fLLp','fL2p']
        for fXX in fXXs:
            for kind in ['onshell','offshell']:
                tab=load('%s/%s-%s-pole-%s.dat'%(self.path2nuctab,fXX,kind,self.group))
                self.D['%s-%s-pole'%(fXX,kind)]=RectBivariateSpline(G,Y,tab)
  
        fXXs=['f22p','fLLp','fL2p','f22n','fLLn','fL2n']
        for fXX in fXXs:
            for kind in ['onshell','offshell']:
                tab=load('%s/%s-%s-cont-%s.dat'%(self.path2nuctab,fXX,kind,self.group))
                self.D['%s-%s-cont'%(fXX,kind)]=RectBivariateSpline(G,Y,tab)
  
    def get_fXX(self,fXX,kind,x,Q2,y):
        gam=(1+4*x**2*self.mN2_/Q2)**0.5
        return self.D['%s-%s'%(fXX,kind)](gam,y)[0,0]
  
    def get_fXX2(self,fXX,kind,x,Q2,y):
        gam=np.array((1+4*x**2*self.mN2_/Q2)**0.5)
        tab = self.tabs[(fXX,kind)]
        return griddata((self.G,self.Y),tab,(gam,y),fill_value=0,method='cubic')

if __name__== "__main__":

    #--choose group 
    helium=HELIUM(path2nuctab='./grids/helium',group='kpsv',precalc=True)
    #helium=HELIUM(path2nuctab='./grids/helium',group='ss',precalc=True)

    #--Generate tables
    helium.set_normalizations()
    helium.gen_tables()



