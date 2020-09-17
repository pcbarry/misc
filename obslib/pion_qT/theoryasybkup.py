#!/usr/bin/env python
import sys,os
import numpy as np
from scipy.integrate import quad, fixed_quad
from scipy.interpolate import interp1d
from scipy.special import gamma
import warnings
warnings.filterwarnings("ignore")
from tools.config import conf
from qcdlib import aux,alphaS,eweak,pdf0
from qcdlib import mellin
from fitlib import parman
import lhapdf
#--local
#import fakepdf
from reader import READER
from tools.parallel import PARALLEL
from tools.tools import save, load, checkdir
from tools.config import conf,load_config,options
#plotting tools
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pylab as py

class PION_QT:

    def __init__(self):
        self.CF=4.0/3.0
        self.alpha = conf['aux'].alfa 
        self.Mp = conf['aux'].M 
        self.eU2 = 4./9
        self.eD2 = 1./9
        self.eq2=np.array([0,self.eU2,self.eU2,self.eD2,self.eD2,self.eD2,self.eD2,self.eU2,self.eU2,self.eU2,self.eU2,self.eD2,self.eD2])
        self.iflav=[21,2,-2,1,-1,3,-3,4,-4,5,-5]
        self.flav = ['g','u','ub','d','db','s','sb','c','cb','b','bb']
        #self.pion_pdf=#pdf0.PDF()
        self.tungsten_pdf=lhapdf.mkPDFs('EPPS16nlo_CT14nlo_W184')

    def integrand1(self,uxa,uy,pT,s,Q,ymax,ipdf):

        #--kinematics
        rs=s**0.5
        Q2=Q**2
        tau=Q2/s
        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2

        #--rapidity
        y = uy*(ymax)
        if np.abs(y)>np.arccosh((1+tau)/(2*np.sqrt((pT**2+Q**2)/s))):
            return 0

        #--pT
        pTmax = rs/2*(1-tau)
        if pT>pTmax: return 0

        #--momentum fractions
        xT=2*pT/rs
        x1=0.5*(xT**2+4*tau)**0.5*np.exp(y)
        x2=0.5*(xT**2+4*tau)**0.5*np.exp(-y)
        xamin=(x1-tau)/(1-x2)
        xa=xamin + uxa*(1-xamin)
        if xa>1: xa = 1
        xb=(xa*x2-tau)/(xa-x1)
        if xb>1: xb = 1

        #--hard factors
        sh=xa*xb*s
        th=-xa*s*x2+Q2
        uh=-xb*s*x1+Q2
        alphas=conf['alphaS'].get_alphaS(mu2)
        # qqb->l,lb,g
        HA=  alphas*self.alpha**2*8./27/(Q2*sh**2)*(th/uh+uh/th+2*Q2*sh/(th*uh))
        # qg ->l,l
        HC= -alphas*self.alpha**2*1./9/(Q2*sh**2)*(th/sh+sh/th+2*Q2*uh/(th*sh))
        HC2 = -alphas*self.alpha**2*1./9/(Q2*sh**2)*(uh/sh+sh/uh+2*Q2*th/(uh*sh))

        #--soft factors
        conf['pdf-pion'].evolve(mu2)
        fA_pion=np.array([conf['pdf-pion'].get_xF(xa,mu2,_,evolve=False)/xa for _ in self.flav])
        fB_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xb,mu2)/xb for _ in self.iflav])

        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)

        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)

        lumA =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[1,3,5,7,9]]*fB_tungsten[[2,4,6,8,10]])
        lumA+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[2,4,6,8,10]]*fB_tungsten[[1,3,5,7,9]])
    
        lumC =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*self.eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[0]*fB_tungsten[[1,2,3,4,5,6,7,8,9,10]])
        lumC2 =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*self.eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[[1,2,3,4,5,6,7,8,9,10]]*fB_tungsten[0])

        dsigma = HA*lumA + HC*lumC + HC2*lumC2
        factor = xa*xb/(xa-x1) * 2 * pT * (2*ymax) * (1-xamin) * 2*Q * 0.3893793656*1e-27
        return dsigma*factor

    def get_xsec1(self,pT,s,Q,ymax):
        """
        returns dsig/dpTdQ  averaged over y
        """
        df2 = lambda uxa,uy: self.integrand1(uxa,uy,pT,s,Q,ymax,0)
        df1 = lambda uy: fixed_quad(lambda uxa: np.vectorize(df2)(uxa,uy),0,1,n=10)[0]
        return fixed_quad(np.vectorize(df1),0,1,n=10)[0]
 
    def integrand2(self,uxa,uQ,pT,xF,s,Qmin,Qmax,ipdf):

        #--kinematics
        Q=Qmin+uQ*(Qmax-Qmin)
        rs=s**0.5
        Q2=Q**2
        tau=Q2/s
        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2
        mT = np.sqrt(Q2+pT**2)

        #--rapidity
        y = np.arcsinh(rs*xF/(2*mT))
        if np.abs(y)>np.arccosh((1+tau)/(2*np.sqrt((pT**2+Q**2)/s))):
            return 0
        xF_jac = (rs)/(2*mT*np.cosh(y))

        #--pT
        pTmax = rs/2*(1-tau)
        if pT>pTmax: return 0

        #--momentum fractions
        xT=2*pT/rs
        x1=0.5*(xT**2+4*tau)**0.5*np.exp(y)
        x2=0.5*(xT**2+4*tau)**0.5*np.exp(-y)
        xamin=(x1-tau)/(1-x2)
        xa=xamin + uxa*(1-xamin)
        if xa>1: xa=1
        xb=(xa*x2-tau)/(xa-x1)
        if xb>1: xb=1

        #--hard factors
        sh=xa*xb*s
        th=-xa*s*x2+Q2
        uh=-xb*s*x1+Q2
        alphas=conf['alphaS'].get_alphaS(mu2)
        # qqb->l,lb,g
        HA=  alphas*self.alpha**2*8./27/(Q2*sh**2)*(th/uh+uh/th+2*Q2*sh/(th*uh))
        # qg ->l,l
        HC= -alphas*self.alpha**2*1./9/(Q2*sh**2)*(th/sh+sh/th+2*Q2*uh/(th*sh))
        HC2 = -alphas*self.alpha**2*1./9/(Q2*sh**2)*(uh/sh+sh/uh+2*Q2*th/(uh*sh))

        #--soft factors
        conf['pdf-pion'].evolve(mu2)
        fA_pion=np.array([conf['pdf-pion'].get_xF(xa,mu2,_,evolve=False)/xa for _ in self.flav])
        fB_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xb,mu2)/xb for _ in self.iflav])
    
        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)

        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)

        lumA =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[1,3,5,7,9]]*fB_tungsten[[2,4,6,8,10]])
        lumA+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[2,4,6,8,10]]*fB_tungsten[[1,3,5,7,9]])
    
        lumC =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*self.eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[0]*fB_tungsten[[1,2,3,4,5,6,7,8,9,10]])
        lumC2 =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*self.eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[[1,2,3,4,5,6,7,8,9,10]]*fB_tungsten[0])

        dsigma = HA*lumA + HC*lumC + HC2*lumC2
        factor = xa*xb/(xa-x1) * 2 * pT * (Qmax-Qmin) * xF_jac * 2*Q * (1-xamin) * 0.3893793656*1e-27

        return dsigma*factor

    def get_xsec2(self,pT,xF,s,Qmin,Qmax):
        df2 = lambda uxa, uQ: self.integrand2(uxa,uQ,pT,xF,s,Qmin,Qmax,0)
        df1 = lambda uQ: fixed_quad(lambda uxa: np.vectorize(df2)(uxa,uQ),0,1,n=10)[0]
        return fixed_quad(np.vectorize(df1),0,1,n=10)[0]

    def asyintRA1(self,Q,pT,s,uy,ipdf):

        rs=s**0.5
        Q2=Q**2
        tau=Q2/s
        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2
        pTmax = ((s**0.5)/2)*(1-(Q**2)/s)
        alfa=1/137. #--should replace with self.alpha
        Mp=0.93891897 #--no need for this I think
        eU2=4./9
        eD2=1./9
        eq2=np.array([0,eU2,eU2,eD2,eD2,eD2,eD2,eU2,eU2,eU2,eU2,eD2,eD2])
        iflav=[21,2,-2,1,-1,3,-3,4,-4,5,-5]

        if pT>pTmax: return 0
        
        ymin = -np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        ymax = np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        
        yuplim = np.arcsinh(s**0.5/(2*np.sqrt(Q2+pT**2)))
        ylowlim = 0
        
        y = ylowlim+uy*(yuplim-ylowlim)
        
        if y>ymax: y = ymax
        if y<ymin: y = ymin
        
        xT=2*pT/rs
        
        xamin = np.exp(y)*Q/rs
        if xamin>1: xamin = 1
        xbmin = np.exp(-y)*Q/rs
        if xbmin>1: xbmin = 1
        
        alphaS=conf['alphaS'].get_alphaS(mu2)
        
        # qqb->l,lb,g
        RA1 = alphaS*alfa**2*8./27/(Q2*s*pT**2)*(2*(-1.5+np.log(Q2/pT**2)))
        
        conf['pdf-pion'].evolve(mu2)
        fA_pion=np.array([conf['pdf-pion'].get_xF(xamin,mu2,_,evolve=False)/xamin for _ in self.flav])
        fB_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xbmin,mu2)/xbmin for _ in self.iflav])
        
        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)

        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)
        
        dsigmaRA1=0
        
        lumA =np.sum(iNf[[1,3,5,7,9]]*eq2[[1,3,5,7,9]]*fA_pion[[1,3,5,7,9]]*fB_tungsten[[2,4,6,8,10]])
        lumA+=np.sum(iNf[[1,3,5,7,9]]*eq2[[1,3,5,7,9]]*fA_pion[[2,4,6,8,10]]*fB_tungsten[[1,3,5,7,9]])
        dsigmaRA1+= RA1*lumA
        
        defaultRA1 = dsigmaRA1 * 2*Q
        
        defaultRA1*=0.3893793656*1e-27
        
        return defaultRA1*2*pT*(ymax-ymin)

    def asyintA(self,Q,uxa,pT,s,uy,ipdf):

        rs=s**0.5
        Q2=Q**2
        tau=Q2/s
        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2
        pTmax = ((s**0.5)/2)*(1-(Q**2)/s)
        alfa=1/137.
        
        if pT>pTmax: return 0
        
        ymin = -np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        ymax = np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        
        yuplim = np.arcsinh(s**0.5/(2*np.sqrt(Q2+pT**2)))
        ylowlim = 0
        
        y = ylowlim+uy*(yuplim-ylowlim)
        
        if y>=ymax: y = ymax
        if y<=ymin: y = ymin
        
        xT=2*pT/rs
        
        x1=0.5*(xT**2+4*tau)**0.5*np.exp(y)
        x2=0.5*(xT**2+4*tau)**0.5*np.exp(-y)
        
        xamin = np.exp(y)*Q/rs
        xbmin = np.exp(-y)*Q/rs
        if xbmin>1: xbmin = 1
        xa=xamin + uxa*(1-xamin)
        if xa>1: xa = 1
        
        alphaS=conf['alphaS'].get_alphaS(mu2)
        
        # qqb->l,lb,g
        RA2 = (1-xamin)*((1+(xamin/xa)**2)/(1-(xamin/xa)))*(alphaS*alfa**2*8./27/(Q2*s*pT**2))
        # qg ->l,l
        RC2 = (1-xamin)*(1/xa)*(1-2*(xamin/xa)+2*(xamin/xa)**2)*(alphaS*alfa**2*1./9/(Q2*s*pT**2))
        
        conf['pdf-pion'].evolve(mu2)
        #fA1_pion=np.array([conf['pdf-pion'].get_xF(xamin/xa,mu2,_,evolve=False)/xa for _ in self.flav])
        fA2_pion=np.array([conf['pdf-pion'].get_xF(xamin,mu2,_,evolve=False)/xa for _ in self.flav])
        fA_pion=np.array([conf['pdf-pion'].get_xF(xa,mu2,_,evolve=False)/xa for _ in self.flav])
        fB_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xbmin,mu2)/xbmin for _ in self.iflav])
        
        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)

        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)
        dsigmaRA2=0
        dsigmaRC2=0
        
        #lumA =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*(fA1_pion[[1,3,5,7,9]]/xa - fA2_pion[[1,3,5,7,9]])*fB_tungsten[[2,4,6,8,10]])
        #lumA+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*(fA1_pion[[2,4,6,8,10]]/xa - fA2_pion[[2,4,6,8,10]])*fB_tungsten[[1,3,5,7,9]])
        lumA =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*(fA_pion[[1,3,5,7,9]]/xa - fA2_pion[[1,3,5,7,9]]*(1 - (0.5*xamin*(2+xamin)-2*np.log(1-xamin))/(1-xamin)))*fB_tungsten[[2,4,6,8,10]])
        lumA+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*(fA_pion[[2,4,6,8,10]]/xa - fA2_pion[[2,4,6,8,10]]*(1 - 0.5*xamin*(2+xamin)-2*np.log(1-xamin))/(1-xamin))*fB_tungsten[[1,3,5,7,9]])
        #lumA =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*(fA_pion[[1,3,5,7,9]]/xa - fA2_pion[[1,3,5,7,9]])*fB_tungsten[[2,4,6,8,10]])
        #lumA+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*(fA_pion[[2,4,6,8,10]]/xa - fA2_pion[[2,4,6,8,10]])*fB_tungsten[[1,3,5,7,9]])
        dsigmaRA2+= RA2*lumA
        
        lumC =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*self.eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[0]*fB_tungsten[[1,2,3,4,5,6,7,8,9,10]])
        #--lumC+=np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[[1,2,3,4,5,6,7,8,9,10]]*fB_tungsten[0])
        dsigmaRC2+= RC2*lumC
        
        defaultRA2 = dsigmaRA2 * 2*Q
        defaultRC2 = dsigmaRC2 * 2*Q
        
        defaultRA2*=0.3893793656*1e-27
        defaultRC2*=0.3893793656*1e-27

        return  defaultRA2*2*pT*(ymax-ymin)+defaultRC2*2*pT*(ymax-ymin)

    def asyintB(self,Q,uxb,pT,s,uy,ipdf):
        rs=s**0.5
        Q2=Q**2
        tau=Q2/s
        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2
        pTmax = ((s**0.5)/2)*(1-(Q**2)/s)
        alfa=1/137.        

        if pT>pTmax: return 0

        ymin = -np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        ymax = np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        
        if (1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s))<1:
              ymin = -np.arccosh(1)
              ymax = np.arccosh(1)
        
        yuplim = np.arcsinh(s**0.5/(2*np.sqrt(Q2+pT**2)))
        ylowlim = 0
        
        y = ylowlim+uy*(yuplim-ylowlim)
        
        if y>=ymax: y = ymax
        if y<=ymin: y = ymin
        
        xT=2*pT/rs
        
        x1=0.5*(xT**2+4*tau)**0.5*np.exp(y)
        x2=0.5*(xT**2+4*tau)**0.5*np.exp(-y)
        
        xbmin = np.exp(-y)*Q/rs
        xamin = np.exp(y)*Q/rs
        if xamin>1: xamin=1
        xb = xbmin + uxb*(1-xbmin)
        if xb>1: xb=1
        
        alphaS=conf['alphaS'].get_alphaS(mu2)
        
        # qqb->l,lb,g
        RA3 = (1-xbmin)*((1+(xbmin/xb)**2)/(1-(xbmin/xb)))*(alphaS*alfa**2*8./27/(Q2*s*pT**2))
        # qg ->l,l
        RC1 = alphaS*alfa**2*1./9/(Q2*s*pT**2)*((1-xbmin)*(1/xb)*(1-2*(xbmin/xb)+2*(xbmin/xb)**2))
        
        conf['pdf-pion'].evolve(mu2)
        fA_pion=np.array([conf['pdf-pion'].get_xF(xamin,mu2,_,evolve=False)/xamin for _ in self.flav])
        fB_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xb,mu2)/xb for _ in self.iflav])
        #fB1_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xbmin/xb,mu2)/xb for _ in self.iflav])
        fB2_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xbmin,mu2)/xb for _ in self.iflav])

        
        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)
        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)
        
        dsigmaRA3=0
        dsigmaRC1=0
        
        #lumA =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[1,3,5,7,9]]*(fB1_tungsten[[2,4,6,8,10]]/xb-fB2_tungsten[[2,4,6,8,10]]))
        #lumA+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[2,4,6,8,10]]*(fB1_tungsten[[1,3,5,7,9]]/xb-fB2_tungsten[[1,3,5,7,9]]))
        lumA =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[1,3,5,7,9]]*(fB_tungsten[[2,4,6,8,10]]/xb-fB2_tungsten[[2,4,6,8,10]]*(1 - (0.5*xbmin*(2+xbmin)-2*np.log(1-xbmin))/(1-xbmin))))
        lumA+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[2,4,6,8,10]]*(fB_tungsten[[1,3,5,7,9]]/xb-fB2_tungsten[[1,3,5,7,9]]*(1 - (0.5*xbmin*(2+xbmin)-2*np.log(1-xbmin))/(1-xbmin))))
        #lumA =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[1,3,5,7,9]]*(fB_tungsten[[2,4,6,8,10]]/xb-fB2_tungsten[[2,4,6,8,10]]))
        #lumA+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion[[2,4,6,8,10]]*(fB_tungsten[[1,3,5,7,9]]/xb-fB2_tungsten[[1,3,5,7,9]]))
        dsigmaRA3+= RA3*lumA

        #--lumC =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[0]*fB_tungsten[[1,2,3,4,5,6,7,8,9,10]])
        lumC=np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*self.eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[[1,2,3,4,5,6,7,8,9,10]]*fB_tungsten[0])
        dsigmaRC1+= RC1*lumC
        
        defaultRA3 = dsigmaRA3 * 2*Q
        defaultRC1 = dsigmaRC1 * 2*Q
        
        defaultRA3*=0.3893793656*1e-27
        defaultRC1*=0.3893793656*1e-27
        
        return defaultRA3*2*pT*(ymax-ymin)+defaultRC1*2*pT*(ymax-ymin)

    def asyint3(self,uQ,Qmin,Qmax,uxF,xFmin,xFmax,pT,s,ipdf):

        Q=Qmin+uQ*(Qmax-Qmin)
        rs=s**0.5
        Q2=Q**2
        tau=Q2/s
        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2
        pTmax = ((s**0.5)/2)*(1-(Q**2)/s)
        alfa=1/137.        

        if pT>pTmax: return 0
        mT = np.sqrt(Q2+pT**2)
        
        xF = xFmin + uxF*(xFmax-xFmin)
        y = np.arcsinh(xF*s**0.5/(2*mT))
        xFJacobian = (s**0.5)/(2*mT*np.cosh(y))
        
        if  y>np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s))):
              y = np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        if y<-np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s))):
              y = -np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        
        xT=2*pT/rs
        
        x1=0.5*(xT**2+4*tau)**0.5*np.exp(y)
        x2=0.5*(xT**2+4*tau)**0.5*np.exp(-y)
        
        xamin = np.exp(y)*Q/rs
        if xamin>1: xamin=1
        xbmin = np.exp(-y)*Q/rs
        if xbmin>1: xbmin=1
        
        alphaS=conf['alphaS'].get_alphaS(mu2)
        
        # qqb->l,lb,g
        RA1 = alphaS*alfa**2*8./27/(Q2*s*pT**2)*(2*(-1.5+np.log(Q2/pT**2)))

        conf['pdf-pion'].evolve(mu2)
        fA_pion=np.array([conf['pdf-pion'].get_xF(xamin,mu2,_,evolve=False)/xamin for _ in self.flav])
        fB_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xbmin,mu2)/xbmin for _ in self.iflav])
        
        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)
        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)
        
        dsigmaRA1=0
        lumA =np.sum(iNf[[1,3,5,7,9]]*eq2[[1,3,5,7,9]]*fA_pion[[1,3,5,7,9]]*fB_tungsten[[2,4,6,8,10]])
        lumA+=np.sum(iNf[[1,3,5,7,9]]*eq2[[1,3,5,7,9]]*fA_pion[[2,4,6,8,10]]*fB_tungsten[[1,3,5,7,9]])
        dsigmaRA1+= RA1*lumA
        
        defaultRA1 = dsigmaRA1 * 2*Q
        defaultRA1*=0.3893793656*1e-27
        
        return defaultRA1*2*pT*(Qmax-Qmin)*xFJacobian

    def asyint3A(self,uQ,Qmin,Qmax,uxa,uxF,xFmin,xFmax,pT,s,ipdf):

        Q=Qmin+uQ*(Qmax-Qmin)
        rs=s**0.5
        Q2=Q**2
        tau=Q2/s
        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2
        pTmax = ((s**0.5)/2)*(1-(Q**2)/s)
        alfa=1/137.        
        Mp=0.93891897
        eU2=4./9
        eD2=1./9
        eq2=np.array([0,eU2,eU2,eD2,eD2,eD2,eD2,eU2,eU2,eU2,eU2,eD2,eD2])
        iflav=[21,2,-2,1,-1,3,-3,4,-4,5,-5]

        if pT>pTmax: return 0
        mT = np.sqrt(Q2+pT**2)
        
        xF = xFmin+ uxF*(xFmax-xFmin)
        y = np.arcsinh(s**0.5*xF/(2*mT))
        xFJacobian = (s**0.5)/(2*mT*np.cosh(y))
        
        if  y>np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s))):
              y = np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        if y<-np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s))):
              y = -np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        
        xT=2*pT/rs
        
        x1=0.5*(xT**2+4*tau)**0.5*np.exp(y)
        x2=0.5*(xT**2+4*tau)**0.5*np.exp(-y)
        
        xamin = np.exp(y)*Q/rs
        xbmin = np.exp(-y)*Q/rs
        if xbmin>1: xbmin=1
        xa=xamin + uxa*(1-xamin)
        if xa>1: xa = 1
        
        alphaS=conf['alphaS'].get_alphaS(mu2)
        
        # qqb->l,lb,g
        RA2 = (1-xamin)*((1+(xamin/xa)**2)/(1-(xamin/xa)))*(alphaS*alfa**2*8./27/(Q2*s*pT**2))
        # qg ->l,l
        RC2 = (1-xamin)*(1/xa)*(1-2*(xamin/xa)+2*(xamin/xa)**2)*(alphaS*alfa**2*1./9/(Q2*s*pT**2))

        conf['pdf-pion'].evolve(mu2)
        fA1_pion=np.array([conf['pdf-pion'].get_xF(xamin/xa,mu2,_,evolve=False)/xa for _ in self.flav])
        fA2_pion=np.array([conf['pdf-pion'].get_xF(xamin,mu2,_,evolve=False)/xa for _ in self.flav])
        fA_pion=np.array([conf['pdf-pion'].get_xF(xa,mu2,_,evolve=False)/xa for _ in self.flav])
        fB_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xbmin,mu2)/xbmin for _ in self.iflav])
        
        
        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)
        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)
        
        dsigmaRA2=0
        dsigmaRC2=0
        
        lumA =np.sum(iNf[[1,3,5,7,9]]*eq2[[1,3,5,7,9]]*(fA1_pion[[1,3,5,7,9]]/xa-fA2_pion[[1,3,5,7,9]])*fB_tungsten[[2,4,6,8,10]])
        lumA+=np.sum(iNf[[1,3,5,7,9]]*eq2[[1,3,5,7,9]]*(fA1_pion[[2,4,6,8,10]]/xa-fA2_pion[[2,4,6,8,10]])*fB_tungsten[[1,3,5,7,9]])
        dsigmaRA2+= RA2*lumA
        
        lumC =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[0]*fB_tungsten[[1,2,3,4,5,6,7,8,9,10]])
        #--lumC+=np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[[1,2,3,4,5,6,7,8,9,10]]*fB_tungsten[0])
        dsigmaRC2+= RC2*lumC
        
        defaultRA2 = dsigmaRA2 * 2*Q
        defaultRC2 = dsigmaRC2 * 2*Q
        
        defaultRA2*=0.3893793656*1e-27
        defaultRC2*=0.3893793656*1e-27
        
        return defaultRA2*2*pT*(Qmax-Qmin)*xFJacobian+defaultRC2*2*pT*(Qmax-Qmin)*xFJacobian

    def asyint3B(self,uQ,Qmin,Qmax,uxb,uxF,xFmin,xFmax,pT,s,ipdf):

        Q=Qmin+uQ*(Qmax-Qmin)
        rs=s**0.5
        Q2=Q**2
        tau=Q2/s
        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2
        pTmax = ((s**0.5)/2)*(1-(Q**2)/s)
        alfa=1/137.        

        if pT>pTmax: return 0
        
        mT = np.sqrt(Q2+pT**2)
        
        xF = xFmin + uxF*(xFmax+xFmin)
        y = np.arcsinh(s**0.5*xF/(2*mT))
        xFJacobian = (s**0.5)/(2*mT*np.cosh(y))
        
        if  y>np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s))):
              y = np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        if y<-np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s))):
              y = -np.arccosh((1+(Q2/s))/(2*np.sqrt((pT**2+Q2)/s)))
        
        xT=2*pT/rs
        
        x1=0.5*(xT**2+4*tau)**0.5*np.exp(y)
        x2=0.5*(xT**2+4*tau)**0.5*np.exp(-y)
        
        xbmin = np.exp(-y)*Q/rs
        xamin = np.exp(y)*Q/rs
        if xamin>1: xamin=1
        xb = xbmin + uxb*(1-xbmin)
        if xb>1: xb = 1
        
        alphaS=conf['alphaS'].get_alphaS(mu2)
        
        # qqb->l,lb,g
        RA3 = (1-xbmin)*((1+(xbmin/xb)**2)/(1-(xbmin/xb)))*(alphaS*alfa**2*8./27/(Q2*s*pT**2))
        # qg ->l,l
        RC1 = alphaS*alfa**2*1./9/(Q2*s*pT**2)*((1-xbmin)*(1/xb)*(1-2*(xbmin/xb)+2*(xbmin/xb)**2))

        conf['pdf-pion'].evolve(mu2)
        fA_pion=np.array([conf['pdf-pion'].get_xF(xamin,mu2,_,evolve=False)/xamin for _ in self.flav])
        fB1_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xbmin/xb,mu2)/xb for _ in self.iflav])
        fB2_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xbmin,mu2)/xb for _ in self.iflav])
        fB_tungsten=np.array([self.tungsten_pdf[ipdf].xfxQ2(_,xb,mu2)/xb for _ in self.iflav])
        
        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)
        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)
        
        dsigmaRA3=0
        dsigmaRC1=0
        
        lumA =np.sum(iNf[[1,3,5,7,9]]*eq2[[1,3,5,7,9]]*fA_pion[[1,3,5,7,9]]*(fB1_tungsten[[2,4,6,8,10]]/xb-fB2_tungsten[[2,4,6,8,10]]))
        lumA+=np.sum(iNf[[1,3,5,7,9]]*eq2[[1,3,5,7,9]]*fA_pion[[2,4,6,8,10]]*(fB1_tungsten[[1,3,5,7,9]]/xb-fB2_tungsten[[1,3,5,7,9]]))
        dsigmaRA3+= RA3*lumA

        #--lumC =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[0]*fB_tungsten[[1,2,3,4,5,6,7,8,9,10]])
        lumC=np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion[[1,2,3,4,5,6,7,8,9,10]]*fB_tungsten[0])
        dsigmaRC1+= RC1*lumC
        
        defaultRA3 = dsigmaRA3 * 2*Q
        defaultRC1 = dsigmaRC1 * 2*Q
        
        defaultRA3*=0.3893793656*1e-27
        defaultRC1*=0.3893793656*1e-27
        
        return defaultRA3*2*pT*(Qmax-Qmin)*xFJacobian+defaultRC1*2*pT*(Qmax-Qmin)*xFJacobian

    def get_asy_Q_term(self,Q,pT,s):
        "This gets the asymptotic term for the dQ dataset"

        #RA= lambda uy,uQ: self.asyintRA1(uQ,Qmin,Qmax,pT,s,uy,0)
        #dRA = lambda uQ: fixed_quad(lambda uy: np.vectorize(RA)(uy,uQ),0,1,n=10)[0]
        RA= lambda uy: self.asyintRA1(Q,pT,s,uy,0)
        asyterm = fixed_quad(np.vectorize(RA),0,1,n=10)[0]

        #A = lambda uxa,uy,uQ: self.asyintA(uQ,Qmin,Qmax,uxa,pT,s,uy,0)
        #da = lambda uy,uQ: fixed_quad(lambda uxa: np.vectorize(A)(uxa,uy,uQ),0,1,n=10)[0]
        #dA = lambda uy: fixed_quad(lambda uQ: np.vectorize(da)(uy,uQ),0,1,n=10)[0]
        A = lambda uxa,uy: self.asyintA(Q,uxa,pT,s,uy,0)
        da = lambda uy: fixed_quad(lambda uxa: np.vectorize(A)(uxa,uy),0,1,n=10)[0]
        asyterm += fixed_quad(np.vectorize(da),0,1,n=10)[0]

        #B = lambda uxb,uy,uQ: self.asyintB(uQ,Qmin,Qmax,uxb,pT,s,uy,0)
        #db = lambda uy,uQ: fixed_quad(lambda uxb: np.vectorize(B)(uxb,uy,uQ),0,1,n=10)[0]
        #dB = lambda uy: fixed_quad(lambda uQ: np.vectorize(db)(uy,uQ),0,1,n=10)[0]
        B = lambda uxb,uy: self.asyintB(Q,uxb,pT,s,uy,0)
        db = lambda uy: fixed_quad(lambda uxb: np.vectorize(B)(uxb,uy),0,1,n=10)[0]
        asyterm += fixed_quad(np.vectorize(db),0,1,n=10)[0]
        
        return asyterm

    def get_asy_xF_term(self,xFmin,xFmax,pT,s,Qmin,Qmax):
        "This gets the asymptotic term for the dxF dataset"

        RA2 = lambda uQ,uxF: self.asyint3(uQ,Qmin,Qmax,uxF,xFmin,xFmax,pT,s,0)
        dRA2 = lambda uxF: fixed_quad(lambda uQ: np.vectorize(RA2)(uQ,uxF),0,1,n=10)[0]
        asyterm2 = fixed_quad(np.vectorize(dRA2),0,1,n=10)[0]

        A2 = lambda uxa,uQ,uxF: self.asyint3A(uQ,Qmin,Qmax,uxa,uxF,xFmin,xFmax,pT,s,0)
        da2 = lambda uQ,uxF: fixed_quad(lambda uxa: np.vectorize(A2)(uxa,uQ,uxF),0,1,n=10)[0]
        dA2 = lambda uxF: fixed_quad(lambda uQ: np.vectorize(da2)(uQ,uxF),0,1,n=10)[0]
        asyterm2 += fixed_quad(np.vectorize(dA2),0,1,n=10)[0]

        B2 = lambda uxb,uQ,uxF: self.asyint3B(uQ,Qmin,Qmax,uxb,uxF,xFmin,xFmax,pT,s,0)
        db2 = lambda uQ,uxF: fixed_quad(lambda uxb: np.vectorize(B2)(uxb,uQ,uxF),0,1,n=10)[0]
        dB2 = lambda uxF: fixed_quad(lambda uQ: np.vectorize(db2)(uQ,uxF),0,1,n=10)[0]
        asyterm2 += fixed_quad(np.vectorize(dB2),0,1,n=10)[0]

        return asyterm2

if __name__=='__main__':

    from qcdlib import alphaS,aux,eweak,mellin,pdfpion0
    from qcdlib import aux,mellin,alphaS,eweak
    
    conf['Q20']   = 1.27**2
    conf['dy_pT_scale'] = 0
    conf['order']='NLO'
    conf['alphaSmode']='backward'
    conf['dglap mode']='truncated'
    conf['scheme']='ZMVFS'
    conf['aux']=aux.AUX()
    conf['alphaS']=alphaS.ALPHAS()
    conf['mellin-pion']=mellin.MELLIN(npts=8,extended=True)
    conf['pdf-pion']=pdfpion0.PDF()
    #pman=parman.PARMAN()
    #pman.set_qpdf_params() 
    conf['datasets']={}
    conf['datasets']['pion_qT']={}
    conf['datasets']['pion_qT']['xlsx']={}
    conf['datasets']['pion_qT']['xlsx'][1001]='pion_qT/expdata/1001.xlsx'
    conf['datasets']['pion_qT']['xlsx'][1002]='pion_qT/expdata/1002.xlsx'
    conf['datasets']['pion_qT']['norm']={}
    conf['datasets']['pion_qT']['filters']=[]
    conf['datasets']['pion_qT']['filters'].append("pT>2.7")
    conf['datasets']['pion_qT']['filters'].append("Q<8.0")
    conf['pion_qT tabs']=READER().load_data_sets('pion_qT')
    
    conf['pion_qT']=PION_QT()

    tabs=conf['pion_qT tabs']

    channels1A_rat = []
    channels1C_rat = []
    channels1C2_rat = []
    channels2A_rat = []
    channels2C_rat = []
    channels2C2_rat = []
    pT_range = []
    pT_range2 = []
    index = 1
    index2 = 1

    for idx in tabs:
        for i in range(len(tabs[idx]['value'])):

            if tabs[idx]['obs'][0] == 'd2sigma/dpTdx':

                fig2 = py.figure(num = 2, figsize = (15,15))
                pT    = tabs[idx]['pT'][i]
                xF    = tabs[idx]['xF'][i]
                s     = tabs[idx]['s'][i]
                Qmin  = tabs[idx]['Qmin'][i]
                Qmax  = tabs[idx]['Qmax'][i]

                if xF < 0.8:
                  if xF == tabs[idx]['xF'][i+1]:
                    thy = conf['pion_qT'].get_xsec2(pT,xF,s,Qmin,Qmax)
                    channels2A_rat.append(np.abs(((conf['pion_qT'].get_xsec2_channels(pT,xF,s,Qmin,Qmax))[0])/thy))
                    channels2C_rat.append(np.abs(((conf['pion_qT'].get_xsec2_channels(pT,xF,s,Qmin,Qmax))[1])/thy))
                    channels2C2_rat.append(np.abs(((conf['pion_qT'].get_xsec2_channels(pT,xF,s,Qmin,Qmax))[2])/thy))
                    pT_range2.append(pT)
                  else:
                    thy = conf['pion_qT'].get_xsec2(pT,xF,s,Qmin,Qmax)
                    channels2A_rat.append(np.abs(((conf['pion_qT'].get_xsec2_channels(pT,xF,s,Qmin,Qmax))[0])/thy))
                    channels2C_rat.append(np.abs(((conf['pion_qT'].get_xsec2_channels(pT,xF,s,Qmin,Qmax))[1])/thy))
                    channels2C2_rat.append(np.abs(((conf['pion_qT'].get_xsec2_channels(pT,xF,s,Qmin,Qmax))[2])/thy))
                    pT_range2.append(pT)

                    #make subplot of given xF bin
                    ax2 = py.subplot(3,3,index2)
                    ax2.plot(pT_range2,channels2A_rat,'r')
                    ax2.plot(pT_range2,channels2C_rat,'g')
                    ax2.plot(pT_range2,channels2C2_rat,'b')
                    ax2.set_title('dsig/dpTdxF versus pT bin '+str(index2))
                    ax2.set_xlabel('pT')
                    ax2.set_ylabel('dsig/dpTdxf')
 
                    channels2A_rat = []
                    channels2C_rat = []
                    channels2C2_rat = []
                    pT_range2 = []
                    index2 += 1

                else:
                  thy = conf['pion_qT'].get_xsec2(pT,xF,s,Qmin,Qmax)
                  channels2A_rat.append(np.abs(((conf['pion_qT'].get_xsec2_channels(pT,xF,s,Qmin,Qmax))[0])/thy))
                  channels2C_rat.append(np.abs(((conf['pion_qT'].get_xsec2_channels(pT,xF,s,Qmin,Qmax))[1])/thy))
                  channels2C2_rat.append(np.abs(((conf['pion_qT'].get_xsec2_channels(pT,xF,s,Qmin,Qmax))[2])/thy))
                  pT_range2.append(pT)  
                  #make subplot of given xF bin                  

                  ax2 = py.subplot(3,3,index2)
                  plot2A = ax2.plot(pT_range2,channels2A_rat,'r',label = 'Channel A')
                  plot2C = ax2.plot(pT_range2,channels2C_rat,'g',label = 'Channel C')
                  plot2C2 = ax2.plot(pT_range2,channels2C2_rat,'b',label = 'Channel C2')
                  ax2.set_title('dsig/dpTdxF versus pT bin '+str(index2))
                  ax2.set_xlabel('pT')
                  ax2.set_ylabel('dsig/dpTdxf')
                  if pT == pT_range[0]:
                    ax2.legend()
                py.tight_layout()
                py.savefig('ratios-xF.png')

            elif tabs[idx]['obs'][0] == 'd2sigma/dpTdm':
                fig = py.figure(num = 1, figsize = (15,15))
                pT   = tabs[idx]['pT'][i]
                s    = tabs[idx]['s'][i]
                Q    = tabs[idx]['Q'][i]
                ymax = tabs[idx]['ymax'][i]

                if Q != 7.2:
                  if Q == tabs[idx]['Q'][i+1]:
                    thy = conf['pion_qT'].get_xsec1(pT,s,Q,ymax)
                    channels1A_rat.append(np.abs(((conf['pion_qT'].get_xsec1_channels(pT,s,Q,ymax))[0])/thy))
                    channels1C_rat.append(np.abs(((conf['pion_qT'].get_xsec1_channels(pT,s,Q,ymax))[1])/thy))
                    channels1C2_rat.append(np.abs(((conf['pion_qT'].get_xsec1_channels(pT,s,Q,ymax))[2])/thy))
                    pT_range.append(pT)
                  else:
                    thy = conf['pion_qT'].get_xsec1(pT,s,Q,ymax)
                    channels1A_rat.append(np.abs(((conf['pion_qT'].get_xsec1_channels(pT,s,Q,ymax))[0])/thy))
                    channels1C_rat.append(np.abs(((conf['pion_qT'].get_xsec1_channels(pT,s,Q,ymax))[1])/thy))
                    channels1C2_rat.append(np.abs(((conf['pion_qT'].get_xsec1_channels(pT,s,Q,ymax))[2])/thy))
                    pT_range.append(pT)

                    #make subplot of given Q bin
                    ax = py.subplot(2,3,index)
                    ax.plot(pT_range,channels1A_rat,'r')
                    ax.plot(pT_range,channels1C_rat,'g')
                    ax.plot(pT_range,channels1C2_rat,'b')
                    ax.set_title('dsig/dpTdQ versus pT bin '+str(index))
                    ax.set_xlabel('pT')
                    ax.set_ylabel('dsig/dpTdQ')

                    channels1A_rat = []
                    channels1C_rat = []
                    channels1C2_rat = []
                    pT_range = []
                    index += 1
        
                else:
                  thy = conf['pion_qT'].get_xsec1(pT,s,Q,ymax)
                  channels1A_rat.append(np.abs(((conf['pion_qT'].get_xsec1_channels(pT,s,Q,ymax))[0])/thy))
                  channels1C_rat.append(np.abs(((conf['pion_qT'].get_xsec1_channels(pT,s,Q,ymax))[1])/thy))
                  channels1C2_rat.append(np.abs(((conf['pion_qT'].get_xsec1_channels(pT,s,Q,ymax))[2])/thy))
                  pT_range.append(pT)

                  #make subplot of given Q bin

                  ax = py.subplot(2,3,index)
                  ax.plot(pT_range,channels1A_rat,'r',label = 'Channel A')
                  ax.plot(pT_range,channels1C_rat,'g',label = 'Channel C')
                  ax.plot(pT_range,channels1C2_rat,'b',label = 'Channel C2')
                  ax.set_title('dsig/dpTdQ versus pT bin '+str(index))
                  ax.set_xlabel('pT')
                  ax.set_ylabel('dsig/dpTdQ')
                  if pT == pT_range[0]:
                    ax.legend()
                py.tight_layout()
                py.savefig('ratios-Q.png')

            print idx,i,thy




