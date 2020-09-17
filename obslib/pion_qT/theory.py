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
        if 'W col' in conf: self.tungsten_pdf=lhapdf.mkPDFs(conf['W col'])
        else: self.tungsten_pdf=lhapdf.mkPDFs('EPPS16nlo_CT14nlo_W184')

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
        y=ymax*(2*uy-1)
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
        if 'W set' in conf: df2 = lambda uxa,uy: self.integrand1(uxa,uy,pT,s,Q,ymax,conf['W set'])
        else: df2 = lambda uxa,uy: self.integrand1(uxa,uy,pT,s,Q,ymax,0)
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
        if 'W set' in conf: df2 = lambda uxa, uQ: self.integrand2(uxa,uQ,pT,xF,s,Qmin,Qmax,conf['W set'])
        else: df2 = lambda uxa, uQ: self.integrand2(uxa,uQ,pT,xF,s,Qmin,Qmax,0)
        df1 = lambda uQ: fixed_quad(lambda uxa: np.vectorize(df2)(uxa,uQ),0,1,n=10)[0]
        return fixed_quad(np.vectorize(df1),0,1,n=10)[0]
    

if __name__=='__main__':

    from qcdlib import alphaS,aux,eweak,mellin,pdfpion0
    from qcdlib import aux,mellin,alphaS,eweak
    
    conf['Q20']   = 1.27**2
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
    #conf['datasets']['pion_qT']['filters'].append("pT>=pTmax")
    conf['pion_qT tabs']=READER().load_data_sets('pion_qT')
    
    conf['pion_qT']=PION_QT()

    tabs=conf['pion_qT tabs']
    for idx in tabs:
        for i in range(len(tabs[idx]['value'])):

            if tabs[idx]['obs'][0] == 'd2sigma/dpTdx':

                pT    = tabs[idx]['pT'][i]
                xF    = tabs[idx]['xF'][i]
                s     = tabs[idx]['s'][i]
                Qmin  = tabs[idx]['Qmin'][i]
                Qmax  = tabs[idx]['Qmax'][i]
                thy = conf['pion_qT'].get_xsec2(pT,xF,s,Qmin,Qmax)

            elif tabs[idx]['obs'][0] == 'd2sigma/dpTdm':

                pT   = tabs[idx]['pT'][i]
                s    = tabs[idx]['s'][i]
                Q    = tabs[idx]['Q'][i]
                ymax = tabs[idx]['ymax'][i]
                thy = conf['pion_qT'].get_xsec1(pT,s,Q,ymax)

            print idx,i,thy





