#!/usr/bin/env python
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma,binom,factorial
from scipy import interpolate
from tools.config import conf
from mpmath import fp

class PITD:

    def __init__(self,pdftype='pion'):
        self.euler=conf['aux'].euler
        self.CF=conf['aux'].CF
        self.aS=conf['alphaS'].get_alphaS
        if pdftype=='pion': self.pdf=conf['pdf-pion-int']
        self.pdftype=pdftype

        #--systematic effects
        self.alpha=-0.5
        self.beta=3
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
        sigma01=interpolate.interp1d(self.NU,self.storage01,fill_value=0,bounds_error=False) #--interpolated function
        sigma02=interpolate.interp1d(self.NU,self.storage02,fill_value=0,bounds_error=False) #--interpolated function
        sig01=sigma01(nu)
        sig02=sigma02(nu)
        #sigma01=self.storage(nu)
        
        #--lattice spacing effects
        ReP=sig01*self.p1 + sig02*self.p2
        ReR=sig01*self.r1 + sig02*self.r2
        latt_spacing_effect = a/np.abs(z)*ReP + a*ReR

        #--higher-twist effect
        ReB=sig01*self.b1 + sig02*self.b2
        higher_twist_effect = z**2*ReB

        #--pion mass effect
        ReM=sig01*self.m1 + sig02*self.b2
        pion_mass_effect = (mpi-conf['aux'].Mpi)*ReM

        #--finite volume effect
        ReF=sig01*self.f1 + sig02*self.f2
        finite_vol_effect = np.exp(-mpi*(L-z))*ReF

        return latt_spacing_effect + higher_twist_effect + pion_mass_effect + finite_vol_effect

    def cn(self,n,z2,mu2):
        CF=self.CF
        euler=self.euler

        #--harmonic number
        def harmonic(ini,fin):
            harmonic=0
            for i in range(ini,fin):
                harmonic+=1./i
            return harmonic

        gamma_n=1.0/((n+1)*(n+2))-0.5-2*harmonic(2,n+1)

        #--polygamma
        def polygamma(ini,fin):
            polygamma=0
            for i in range(ini,fin):
                polygamma+=1./i**2
            return polygamma

        d_n=2*(harmonic(1,n)**2+polygamma(1,n)+0.5-1.0/((n+1)*(n+2)))
        aS=self.aS(mu2)
        cn=1.-aS/2./np.pi*CF*(gamma_n*np.log(z2*mu2*np.exp(2*euler+1)/4.0)+d_n)
        return cn

    def ReITD(self,nu,z2,mpi,L,a):
        z2=25*z2 #--units from fm^2 to GeV^-2
        a=5*a #--units from fm to GeV^-1
        mpi=1e-3*mpi #--units from MeV to GeV
        L=L*a #--from units of a to GeV^-1

        if conf['LCS scale']==0: mu2=2.0**2
        elif conf['LCS scale']==1: mu2=4.0**2
        elif conf['LCS scale']==2: mu2=conf['Q20']
        elif type(conf['LCS scale'])!=int: mu2=conf['LCS scale']**2*conf['Q20']

        N=conf['imell'].N
        M=1
        for n in range(1,len(N)):
            if n%2==1: continue
            self.pdf.evolve(mu2)
            if self.pdftype=='pion': val=np.array(self.pdf.storage[mu2]['ub'])-np.array(self.pdf.storage[mu2]['u'])
            M+=(1j*nu)**n/int(fp.factorial(n))*self.cn(n+1,z2,mu2)*val[n] #--nth index in pdf storage is N=n+1

        #if 'FV nu exp - pITD' in conf: power=conf['FV nu exp - pITD']
        #else: power=0
        #correction=self.FV_mod*np.exp(-mpi*(L-z2**0.5))*nu**2 #--finite volume correction
        correction=self.get_correction(z2**0.5,a,mpi,L,nu)

        M+=correction

        return np.real(M)

    #--get, set state

    def get_state(self):
        return (self.storage01,self.storage02,self.b1,self.p1,self.r1,self.m1,self.f1,self.b2,self.p2,self.r2,self.m2,self.f2)

    def set_state(self,state):
        self.storage01,self.storage02,self.b1,self.p1,self.r1,self.m1,self.f1,self.b2,self.p2,self.r2,self.m2,self.f2=state[:]

    #--parallelization

    def setup_parallel(self):
        self.tasks=[]
        for k in conf['pITD tabs']:
            tab=conf['pITD tabs'][k]
            for i in range(len(tab['obs'])):
                task={}
                task['reaction']='pITD'
                task['obs']=tab['obs'][i]
                task['nu']=tab['nu'][i]
                task['z2']=tab['z'][i]**2
                task['mpi']=tab['mpi'][i]
                task['L']=tab['L'][i]
                task['a']=tab['a'][i]
                self.tasks.append(task)
        self.thyval={}
        for k in conf['pITD tabs']:
            tab=conf['pITD tabs'][k]
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
        task['value']=self.ReITD(nu,z2,mpi,L,a)
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

    def get_ReITD_values(self,nu,z2,mpi,L,a):
        key=(nu,z2,mpi,L,a)
        return self.thyval[key]

if __name__=="__main__":

    from qcdlib import mellin,alphaS,aux,pdfpion

    conf['aux']=aux.AUX()
    conf['order']='NLO'
    conf['Q20']=1.27**2
    conf['dglap mode']='truncated'
    conf['alphaS']=alphaS.ALPHAS()

    imell=mellin.IMELLIN()
    imell.N=np.array([])
    for i in range(1,51):
        iN=np.array([i],dtype=complex)
        imell.N=np.append(imell.N,iN)

    conf['imell']=imell
    conf['pdf-pion-int']=pdfpion.PDF(mellin=conf['imell'])
    conf['LCS scale']=0

    pitd=PITD()
    nu=0.1
    z2=1.0/1.27**2
    mpi=435
    L=32
    a=0.127
    print(pitd.ReITD(nu,z2,mpi,L,a))
