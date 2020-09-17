#!/usr/bin/env python
from scipy.integrate import fixed_quad
from scipy.special import gamma
import numpy as np
from tools.config import conf

class DY:
  
    def __init__(self):
        self.aux=conf['aux']
        if conf['order']=='LO':  self.iord=0
        if conf['order']=='NLO': self.iord=1
        self.mellin=conf['mellin']
        self.eweak=conf['eweak']
        self.CF=self.aux.CF
        self.TF=self.aux.TF
        self.eU2=4./9.0
        self.eD2=1./9.0
  
    def beta(self,a,b):
        return gamma(a)*gamma(b)/gamma(a+b)
  
    def quad(self,f,xmin,xmax):
        #return _quad(f,xmin,xmax)[0]
        return fixed_quad(f,xmin,xmax,n=5)[0]
  
    def get_x1(self,z,y):
        out=self.tau/z*(1-(1-y)*(1-z))/(1-y*(1-z))
        return out**0.5*np.exp(self.Y)
   
    def get_x2(self,z,y):
        out=self.tau/z*(1-y*(1-z))/(1-(1-y)*(1-z))
        return out**0.5*np.exp(-self.Y)
  
    def _get_lum(self,fA,fB,channel): 
  
        if  channel=='qA,qbB':
  
            out= self.eU2*(fA[1]*fB[2])\
                +self.eD2*(fA[3]*fB[4])\
                +self.eD2*(fA[5]*fB[6])
            if self.Nf==4: out+=self.eU2*(fA[7]*fB[8])
            if self.Nf==5: out+=self.eU2*(fA[7]*fB[8])+self.eD2*(fA[9]*fB[10])
  
        elif channel=='qbA,qB':
  
            out= self.eU2*(fB[1]*fA[2])\
                +self.eD2*(fB[3]*fA[4])\
                +self.eD2*(fB[5]*fA[6])
            if self.Nf==4: out+=self.eU2*(fB[7]*fA[8])
            if self.Nf==5: out+=self.eU2*(fB[7]*fA[8])+self.eD2*(fB[9]*fA[10])
  
        elif channel=='qA,gB':
  
            out= self.eU2*(fA[1]+fA[2])*fB[0]\
                +self.eD2*(fA[3]+fA[4])*fB[0]\
                +self.eD2*(fA[5]+fA[6])*fB[0]
            if self.Nf==4: out+=self.eU2*(fA[7]+fA[8])*fB[0] 
            if self.Nf==5: out+=self.eU2*(fA[7]+fA[8])*fB[0]+self.eD2*(fA[9]+fA[10])*fB[0]
  
        elif channel=='gA,qB':
  
            out= self.eU2*(fB[1]+fB[2])*fA[0]\
                +self.eD2*(fB[3]+fB[4])*fA[0]\
                +self.eD2*(fB[5]+fB[6])*fA[0]
            if self.Nf==4: out+=self.eU2*(fB[7]+fB[8])*fA[0]
            if self.Nf==5: out+=self.eU2*(fB[7]+fB[8])*fA[0]+self.eD2*(fB[9]+fB[10])*fA[0]
  
        return out 
  
    def _get_lum_flav(self,fA,fB,flav): 
  
        if   flav=='g' : return fA*fB[0]
        elif flav=='u' : return fA*fB[1]
        elif flav=='ub': return fA*fB[2]
        elif flav=='d' : return fA*fB[3]
        elif flav=='db': return fA*fB[4]
        elif flav=='s' : return fA*fB[5]
        elif flav=='sb': return fA*fB[6]
        elif flav=='c' : return fA*fB[7]
        elif flav=='cb': return fA*fB[8]
        elif flav=='b' : return fA*fB[9]
        elif flav=='bb': return fA*fB[10]
  
    def get_lum(self,z,y,channel): 
        x1=self.get_x1(z,y)
        x2=self.get_x2(z,y)
        if x1>1 or x2>1: return 0  
  
        if  self.ilum=='normal':
            fA=conf['pdfA'].get_pdfs(x1,self.muF2)
            fB=conf['pdfB'].get_pdfs(x2,self.muF2)
            return self._get_lum(fA,fB,channel)
  
        elif self.ilum=='mell-real':
  
            return np.real(x1**(-self.n) * x2**(-self.m))
  
        elif self.ilum=='mell-imag':
  
            return np.imag(x1**(-self.n) * x2**(-self.m))
  
        elif self.ilum=='hybrid-real':
  
            fA=np.real(x1**(-self.n))
            fB=conf['pdfB'].get_pdfs(x2,self.muF2)
            return self._get_lum_flav(fA,fB,self.flav)
  
        elif self.ilum=='hybrid-imag':
  
            fA=np.imag(x1**(-self.n))
            fB=conf['pdfB'].get_pdfs(x2,self.muF2)
            return self._get_lum_flav(fA,fB,self.flav)
  
    def get_jac(self,z,y):
        return 1.0/(1-y*(1-z))/(1-(1-y)*(1-z))
  
    #--qqb channel
  
    #--part 0
  
    def Cqqb0(self,f):
        return (1+self.iord*self.CF*self.aS/np.pi*(3/2.*np.log(self.Q2/self.muF2)+2*np.pi**2/3-6))*0.5*(f(0)+f(1))
    
    def xsecCqqb0AB(self):
        return  self.Cqqb0(lambda y: self.get_lum(1,y,'qA,qbB')*self.get_jac(1,y))
  
    def xsecCqqb0BA(self):
        return  self.Cqqb0(lambda y: self.get_lum(1,y,'qbA,qB')*self.get_jac(1,y))
  
    #--part 1
  
    def Cqqb1(self,z,f):
        plus =np.log(self.Q2/self.muF2)*(1/(1-z)*((1+z*z)*f(z)-2*f(1)) + 2*f(1)*np.log(1-self.tau)/(1-self.tau))
        plus+=2*np.log(1-z)/(1-z)*((1+z*z)*f(z)-2*f(1)) + 2*f(1)*np.log(1-self.tau)**2/(1-self.tau)
        plus+=-1/(1-z)*(1+z*z)*f(z)*np.log(z) + 2*f(1)*(1-np.pi**2/6)/(1-self.tau)
        return 0.5*self.iord*self.CF*self.aS/np.pi*(plus+(1-z)*f(z))
  
    def dxsecCqqb1AB(self,z):
        return  self.Cqqb1(z,lambda z:self.get_lum(z,0,'qA,qbB')*self.get_jac(z,0)+self.get_lum(z,1,'qA,qbB')*self.get_jac(z,1))
  
    def dxsecCqqb1BA(self,z):
        return  self.Cqqb1(z,lambda z:self.get_lum(z,0,'qbA,qB')*self.get_jac(z,0)+self.get_lum(z,1,'qbA,qB')*self.get_jac(z,1))
  
    def xsecCqqb1AB(self):
        return self.quad(np.vectorize(self.dxsecCqqb1AB),self.tau,1)
  
    def xsecCqqb1BA(self):
        return self.quad(np.vectorize(self.dxsecCqqb1BA),self.tau,1)
  
    #--part 2
  
    def aux1(self,z,y):
        return 0.5*(1+(1-z)**2/z*y*(1-y))
  
    def Cqqb2(self,z,y,f):
        plus10=1.0/y*(f(z,y)*self.aux1(z,y)-f(z,0)*self.aux1(z,0))
        plus00=1.0/y*(f(1,y)*self.aux1(1,y)-f(1,0)*self.aux1(1,0))
        plus11=1.0/(1-y)*(f(z,y)*self.aux1(z,y)-f(z,1)*self.aux1(z,1))
        plus01=1.0/(1-y)*(f(1,y)*self.aux1(1,y)-f(1,1)*self.aux1(1,1))
        return self.iord*self.CF*self.aS/np.pi\
              *(  1/(1-z)*((1+z*z)*(plus10+plus11)-2*(plus00+plus01))\
                 +2*np.log(1-self.tau)*(plus00+plus01)/(1-self.tau)\
                 -2*(1-z)*f(z,y)* self.aux1(z,y))
  
    def dxsecCqqb2AB(self,z,y):
        return self.Cqqb2(z,y,lambda z,y: self.get_lum(z,y,'qA,qbB')*self.get_jac(z,y))
  
    def dxsecCqqb2BA(self,z,y):
        return self.Cqqb2(z,y,lambda z,y: self.get_lum(z,y,'qbA,qB')*self.get_jac(z,y))
  
    def xsecCqqb2AB(self):
        func=lambda z: self.quad(np.vectorize(lambda y: self.dxsecCqqb2AB(z,y)),0,1)
        return self.quad(np.vectorize(func),self.tau,1)
  
    def xsecCqqb2BA(self):
        func=lambda z: self.quad(np.vectorize(lambda y: self.dxsecCqqb2BA(z,y)),0,1)
        return self.quad(np.vectorize(func),self.tau,1)
  
    #--combined
  
    def xsecCqqbAB(self):
        return self.xsecCqqb0AB()+self.xsecCqqb1AB()+self.xsecCqqb2AB()
  
    def xsecCqqbBA(self):
        return self.xsecCqqb0BA()+self.xsecCqqb1BA()+self.xsecCqqb2BA()
  
    def xsecCqqb(self):
        return self.xsecCqqbAB() + self.xsecCqqbBA()
  
    #--qg channel
  
    #--part 1
  
    def Cqg1(self,z): 
        return self.iord*self.TF*self.aS/2/np.pi\
            *( (z**2+(1-z)**2) *(np.log(self.Q2/self.muF2)+np.log((1-z)**2/z))+2*z*(1-z))
  
    def dxsecCqg1AB(self ,z): 
        return  self.Cqg1(z)*self.get_lum(z,0,'qA,gB')*self.get_jac(z,0)
  
    def xsecCqg1AB(self):
        return self.quad(np.vectorize(self.dxsecCqg1AB),self.tau,1)
  
    def dxsecCqg1BA(self ,z): 
        return  self.Cqg1(z)*self.get_lum(z,1,'gA,qB')*self.get_jac(z,1)
  
    def xsecCqg1BA(self):
        return self.quad(np.vectorize(self.dxsecCqg1BA),self.tau,1)
  
    #--part 2
  
    def aux2(self,z,y):
        return 1+(1-z)**2/z*y*(1-y)
  
    def Cqg2AB(self,z,y,f):
        plus=1.0/y*(f(z,y)*self.aux2(z,y)-f(z,0)*self.aux2(z,0))
        return self.iord*self.TF*self.aS/2/np.pi\
              *( (z*z+(1-z)**2)*plus+(2*z*(1-z)+(1-z)**2*y)*f(z,y)*self.aux2(z,y))
  
    def Cqg2BA(self,z,y,f):
        plus=1.0/(1-y)*(f(z,y)*self.aux2(z,1-y)-f(z,1)*self.aux2(z,0))
        return self.iord*self.TF*self.aS/2/np.pi\
            *( (z*z+(1-z)**2)*plus+(2*z*(1-z)+(1-z)**2*(1-y))*f(z,y)*self.aux2(z,1-y))
  
    def dxsecCqg2AB(self,z,y):
        return self.Cqg2AB(z,y,lambda z,y: self.get_lum(z,y,'qA,gB')*self.get_jac(z,y))\
  
    def dxsecCqg2BA(self,z,y):
        return self.Cqg2BA(z,y,lambda z,y: self.get_lum(z,y,'gA,qB')*self.get_jac(z,y))
  
    def xsecCqg2AB(self):
        func=lambda z: self.quad(np.vectorize(lambda y: self.dxsecCqg2AB(z,y)),0,1)
        return self.quad(np.vectorize(func),self.tau,1)
  
    def xsecCqg2BA(self):
        func=lambda z: self.quad(np.vectorize(lambda y: self.dxsecCqg2BA(z,y)),0,1)
        return self.quad(np.vectorize(func),self.tau,1)
  
    #--combined
  
    def xsecCqgAB(self):
        return self.xsecCqg1AB()+self.xsecCqg2AB()
  
    def xsecCqgBA(self):
        return self.xsecCqg1BA()+self.xsecCqg2BA()
  
    def xsecCqg(self):
        return self.xsecCqgAB()+self.xsecCqgBA()
  
    #--combined all
  
    def get_xsec(self,Q2,S,Y,muF2,ilum='normal',part='full'):
        self.ilum=ilum
        self.Q2=Q2
        self.muF2=muF2
        self.tau=Q2/S
        self.Y=Y
        self.aS=conf['alphaS'].get_alphaS(self.muF2)
        self.Nf=conf['alphaS'].get_Nf(muF2)
        aEM=self.eweak.get_alpha(Q2)
        factor=(4*np.pi*aEM**2)/(9*Q2*S)
  
        if   part=='qA,qbB': xsec=self.xsecCqqbAB()
        elif part=='qbA,qB': xsec=self.xsecCqqbBA()
        elif part=='q,qb':   xsec=self.xsecCqqb()
        elif part=='qA,gB':  xsec=self.xsecCqgAB()
        elif part=='gA,qB':  xsec=self.xsecCqgBA()
        elif part=='q,g':    xsec=self.xsecCqg()
        elif part=='full':   xsec=self.xsecCqqb()+self.xsecCqg()
  
        return factor*xsec
  
    def _get_xsec(self,Q2,S,Y,muF2,ixsec=2,units='pb',ilum='xspace'):
        xsec=self.get_xsec(Q2,S,Y,muF2)
        #xsec*=self.norm*self.aEM**2/Q2/S
        if   units=='mb': xsec*=0.3894
        elif units=='pb': xsec*=0.3894*1e9
        if   ixsec==1: return xsec              # dsig/dM2dY
        elif ixsec==2: return xsec*2*Q2**0.5    # dsig/dMdY
  
if __name__=='__main__':

    import os
    from qcdlib import aux,mellin,alphaS,eweak,pdf0
    from fakepdf import FAKEPDF

    conf['Q20']   = 1.0
    conf['alphaSmode']='backward'
    conf['order']='NLO'
    conf['aux']=aux.AUX()
    conf['mellin']=mellin.MELLIN(npts=4)
    conf['alphaS']=alphaS.ALPHAS()
    conf['eweak']=eweak.EWEAK()
    conf['pdfA']=FAKEPDF()
    conf['pdfB']=FAKEPDF()

    #--simple tests using fake PDFs

    Q2=10.0
    S=32.0**2
    Y=0.0
    muF2=Q2
    dy=DY()
    print dy.get_xsec(Q2,S,Y,muF2,ilum='normal',part='full')












