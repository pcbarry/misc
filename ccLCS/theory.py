#!/usr/bin/env python
import numpy as np
from tools.config import conf
from mpmath import fp
from scipy.integrate import fixed_quad,quad
from scipy.special import gamma,binom,factorial
from scipy import interpolate

class CCLCS:

    def __init__(self,pdftype='pion'):
        self.euler=conf['aux'].euler
        self.CF=conf['aux'].CF
        self.aS=conf['alphaS'].get_alphaS
        if pdftype=='pion': self.pdf=conf['pdf-pion']
        self.pdftype=pdftype
        #--systematic effects
        self.alpha=-0.5
        self.beta=3
        self.b0=0.0
        self.p0=0.0
        self.r0=0.0
        self.m0=0.0
        self.f0=0.0
        self.b1=0.0
        self.p1=0.0
        self.r1=0.0
        self.m1=0.0
        self.f1=0.0
        self.b2=0.0
        self.p2=0.0
        self.r2=0.0
        self.m2=0.0
        self.f2=0.0
        self.setup()

    def setup(self):
        self.NU=np.linspace(0.19,4.75,40)
        self.storage00=np.array([self.get_sigma0n(0,_) for _ in self.NU])
        self.storage01=np.array([self.get_sigma0n(1,_) for _ in self.NU])
        self.storage02=np.array([self.get_sigma0n(2,_) for _ in self.NU])

    def jacobipoly(self,n,x):
        alpha,beta=self.alpha,self.beta
        jacobi=0
        for j in range(n+1):
            omega=binom(n,j)*(-1)**j/factorial(n)*gamma(alpha+n+1)*gamma(alpha+beta+n+j+1)/gamma(alpha+beta+n+1)/gamma(alpha+j+1)
            jacobi+=omega*x**j
        return jacobi

    def get_sigma0n(self,n,nu):
        #--for real cross section
        alpha,beta=self.alpha,self.beta
        integrand=lambda x: np.cos(nu*x)*x**alpha*(1-x)**beta*self.jacobipoly(n,x)
        return quad(integrand,0,1)[0]

    def get_correction(self,z,a,mpi,L,nu):
        #--real corrections
        sigma00=interpolate.interp1d(self.NU,self.storage00,fill_value=0,bounds_error=False) #--interpolated sigma00 function
        sigma01=interpolate.interp1d(self.NU,self.storage01,fill_value=0,bounds_error=False) #--interpolated sigma01 function
        sigma02=interpolate.interp1d(self.NU,self.storage02,fill_value=0,bounds_error=False) #--interpolated sigma02 function
        sig00=sigma00(nu)
        sig01=sigma01(nu)
        sig02=sigma02(nu)

        #--lattice spacing effects
        ReP=sig00*self.p0 + sig01*self.p1 + sig02*self.p2
        ReR=sig00*self.r0 + sig01*self.r1 + sig02*self.r2
        latt_spacing_effect = a/np.abs(z)*ReP + a*ReR

        #--higher-twist effect
        ReB=sig00*self.b0 + sig01*self.b1 + sig02*self.b2
        higher_twist_effect = z**2*ReB

        #--pion mass effect
        ReM=sig00*self.m0 + sig01*self.m1 + sig02*self.m2
        pion_mass_effect = (mpi-conf['aux'].Mpi)*ReM

        #--finite volume effect
        ReF=sig00*self.f0 + sig01*self.f1 + sig02*self.f2
        finite_vol_effect = np.exp(-mpi*(L-z))*ReF

        return latt_spacing_effect + higher_twist_effect + pion_mass_effect + finite_vol_effect

    def get_y_integral(self,x,nu,mu2,z2):
        aS=self.aS(mu2)
        TLO=lambda y: np.cos(x*nu)
        T10=lambda y: aS*self.CF/2/np.pi*(0.5*np.cos(y*nu)-(2*np.log(1-y)/(1-y)-(y**2-3*y+1)/(1-y))*(np.cos(y*x*nu)-np.cos(x*nu)))
        T11=lambda y: -aS*self.CF/2/np.pi*np.log(z2*mu2/4.0*np.exp(2*self.euler))*(1+y**2)/(1-y)*(np.cos(y*x*nu)-np.cos(x*nu))
        return fixed_quad(lambda y: TLO(y)+T10(y)+T11(y),0,1,n=100)[0]
        #return fixed_quad(lambda y: TLO(y),0,1,n=100)[0]

    def get_x_integral(self,nu,mu2,z2):
        self.pdf.evolve(mu2)
        if self.pdftype=='pion': val=lambda x: self.pdf.get_xF(x,mu2,'ub')/x - self.pdf.get_xF(x,mu2,'u')/x #--valence in terms of pi^-
        integrand=lambda x: val(x)*self.get_y_integral(x,nu,mu2,z2)
        return fixed_quad(np.vectorize(integrand),0,1,n=100)[0]

    def current_current(self,nu,z2,mpi,L,a):
        z2=25*z2 #--units from fm^2 to GeV^-2
        a=5*a #--units from fm to GeV^-1
        mpi=1e-3*mpi #--units from MeV to GeV
        LGeVinv=L*a #--from units of lattice spacing to GeV^-1

        if conf['LCS scale']==0: mu2=2.0**2
        elif conf['LCS scale']==1: mu2=4.0**2
        elif conf['LCS scale']==2: mu2=conf['Q20']
        elif type(conf['LCS scale'])!=int: mu2=conf['LCS scale']**2*conf['Q20']

        T1=1.0/np.pi**2 * self.get_x_integral(nu,mu2,z2)

        ##--finite volume effect
        #if 'FV nu exp - ccLCS' in conf: power=conf['FV nu exp - ccLCS']
        #else: power=0
        #FVcorrection=self.FV_mod*np.exp(-mpi*(L-z2**0.5))*nu**power

        ##--pion mass effect
        #mpicorrection=self.mpi_mod*(mpi-conf['aux'].Mpi)

        ##--lattice spacing effect
        #acorrection=self.a_mod*a

        correction=self.get_correction(z2**0.5,a,mpi,L,nu)
        T1+=correction

        return np.real(T1)

    #--get, set state

    def get_state(self):
        return (self.storage00,self.storage01,self.storage02,self.b0,self.p0,self.r0,self.m0,self.f0,self.b1,self.p1,self.r1,self.m1,self.f1,self.b2,self.p2,self.r2,self.m2,self.f2)

    def set_state(self,state):
        self.storage00,self.storage01,self.storage02,self.b0,self.p0,self.r0,self.m0,self.f0,self.b1,self.p1,self.r1,self.m1,self.f1,self.b2,self.p2,self.r2,self.m2,self.f2=state[:]

    #--parallelization

    def setup_parallel(self):
        self.tasks=[]
        for k in conf['ccLCS tabs']:
            tab=conf['ccLCS tabs'][k]
            for i in range(len(tab['obs'])):
                task={}
                task['reaction']='ccLCS'
                task['obs']=tab['obs'][i]
                task['nu']=tab['nu'][i]
                task['z2']=tab['z'][i]**2
                task['mpi']=tab['mpi'][i]
                task['L']=tab['L'][i]
                task['a']=tab['a'][i]
                self.tasks.append(task)
        self.thyval={}
        for k in conf['ccLCS tabs']:
            tab=conf['ccLCS tabs'][k]
            for i in range(len(tab['obs'])):
                obs=tab['obs'][i]
                nu =tab['nu'][i]
                z  =tab['z'][i]
                z2=z**2
                mpi=tab['mpi'][i]
                L  =tab['L'][i]
                a  =tab['a'][i]
                key=(nu,z2,mpi,L,a)
                self.thyval[key]=0.0

    def process_request(self,task):
        nu=task['nu']
        z2=task['z2']
        mpi=task['mpi']
        L=task['L']
        a=task['a']
        task['value']=self.current_current(nu,z2,mpi,L,a)
        return task

    def update(self,task):
        nu=task['nu']
        z2=task['z2']
        mpi=task['mpi']
        L=task['L']
        a=task['a']
        key=(nu,z2,mpi,L,a)
        self.thyval[key]=task['value']

    def get_tasks(self):
        return self.tasks

    #--master functions to be called by residuals

    def get_ccLCS(self,nu,z2,mpi,L,a):
        key=(nu,z2,mpi,L,a)
        return self.thyval[key]

if __name__=="__main__":

    from qcdlib import mellin,alphaS,aux,pdfpion

    conf['aux']=aux.AUX()
    conf['order']='NLO'
    conf['Q20']=1.27**2
    conf['dglap mode']='truncated'
    conf['alphaS']=alphaS.ALPHAS()
    conf['mellin-pion']=mellin.MELLIN(npts=8,extended=True)

    conf['pdf-pion']=pdfpion.PDF()

    conf['LCS scale']=0
    ccLCS=CCLCS()
    nu=0.1
    z2=1.0/1.27**2
    print(ccLCS.current_current(nu,z2,0.415,32,0.127))
