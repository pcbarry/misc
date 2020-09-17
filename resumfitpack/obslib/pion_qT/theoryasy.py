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
        self.tungsten_pdf=lhapdf.mkPDF('EPPS16nlo_CT14nlo_W184',0)

    def Rqqb1zazb(self,Q,pT,s,y,za,zb):

        Q2=Q**2
        tau=Q2/s

        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2

        #--rapidity
        if np.abs(y)>np.arccosh((1+tau)/(2*np.sqrt((pT**2+Q**2)/s))):
            return 0

        #--pT
        pTmax = s**0.5/2*(1-tau)
        if pT>pTmax: return 0

        alphaS=conf['alphaS'].get_alphaS(mu2)

        xa=tau**0.5*np.exp(y)
        xb=tau**0.5*np.exp(-y)

        prefactor = 8*self.alpha**2*alphaS/(27.0*Q2*s*pT**2)

        fA_pion0=np.array([conf['pdf-pion'].get_xF(xa,mu2,_)/xa for _ in self.flav])
        fB_tungsten0=np.array([self.tungsten_pdf.xfxQ2(_,xb,mu2)/xb for _ in self.iflav])

        fA_pionZ=np.array([conf['pdf-pion'].get_xF(xa/za,mu2,_)/(xa/za) for _ in self.flav])
        fB_tungstenZ=np.array([self.tungsten_pdf.xfxQ2(_,xb/zb,mu2)/(xb/zb) for _ in self.iflav])
        

        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)

        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)

        lum0 =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion0[[1,3,5,7,9]]*fB_tungsten0[[2,4,6,8,10]])
        lum0+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion0[[2,4,6,8,10]]*fB_tungsten0[[1,3,5,7,9]])
        lumZA =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pionZ[[1,3,5,7,9]]*fB_tungsten0[[2,4,6,8,10]])
        lumZA+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pionZ[[2,4,6,8,10]]*fB_tungsten0[[1,3,5,7,9]])
        lumZB =np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion0[[1,3,5,7,9]]*fB_tungstenZ[[2,4,6,8,10]])
        lumZB+=np.sum(iNf[[1,3,5,7,9]]*self.eq2[[1,3,5,7,9]]*fA_pion0[[2,4,6,8,10]]*fB_tungstenZ[[1,3,5,7,9]])

        Rqqb1=0
        #--line 1
        Rqqb1+=2*(np.log(Q**2/pT**2)-3.0/2.0)*lum0/(1-xa)/(1-xb)

        #--line 2
        Rqqb1+=(1+zb**2)/(1-zb)*(lumZB/zb-lum0)/(1-xa)
        Rqqb1+=lum0/(1-xa)/(1-xb)*(1.0/2.0*xb*(2+xb)+2*np.log(1-xb))

        #--line 3
        Rqqb1+=(1+za**2)/(1-za)*(lumZA/za-lum0)/(1-xb)
        Rqqb1+=lum0/(1-xa)/(1-xb)*(1.0/2.0*xa*(2+xa)+2*np.log(1-xa))

        return Rqqb1*prefactor

    def Rqg1zazb(self,Q,pT,s,y,za,zb):

        Q2=Q**2
        tau=Q2/s

        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2

        #--rapidity
        if np.abs(y)>np.arccosh((1+tau)/(2*np.sqrt((pT**2+Q**2)/s))):
            return 0

        #--pT
        pTmax = s**0.5/2*(1-tau)
        if pT>pTmax: return 0

        alphaS=conf['alphaS'].get_alphaS(mu2)

        xa=tau**0.5*np.exp(y)
        xb=tau**0.5*np.exp(-y)

        prefactor = self.alpha**2*alphaS/(9.0*Q2*s*pT**2)

        fA_pion0=np.array([conf['pdf-pion'].get_xF(xa,mu2,_)/xa for _ in self.flav])
        fB_tungstenZ=np.array([self.tungsten_pdf.xfxQ2(_,xb/zb,mu2)/(xb/zb) for _ in self.iflav])
        

        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)

        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)

        lumZB  =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*self.eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pion0[[1,2,3,4,5,6,7,8,9,10]]*fB_tungstenZ[[0]])

        Rqg1=(zb**2+(1-zb)**2)/zb/pT**2*lumZB/(1-xa)

        return Rqg1*prefactor

    def Rgq1zazb(self,Q,pT,s,y,za,zb):

        Q2=Q**2
        tau=Q2/s

        if conf['dy_pT_scale'] == 0: mu=Q
        elif conf['dy_pT_scale'] == 1: mu=pT/2
        elif conf['dy_pT_scale'] == 2: mu=pT
        elif conf['dy_pT_scale'] == 3: mu=2*pT
        mu2=mu**2

        #--rapidity
        if np.abs(y)>np.arccosh((1+tau)/(2*np.sqrt((pT**2+Q**2)/s))):
            return 0

        #--pT
        pTmax = s**0.5/2*(1-tau)
        if pT>pTmax: return 0

        alphaS=conf['alphaS'].get_alphaS(mu2)

        xa=tau**0.5*np.exp(y)
        xb=tau**0.5*np.exp(-y)

        prefactor = self.alpha**2*alphaS/(9.0*Q2*s*pT**2)

        fA_pionZ=np.array([conf['pdf-pion'].get_xF(xa/za,mu2,_)/(xa/za) for _ in self.flav])
        fB_tungsten0=np.array([self.tungsten_pdf.xfxQ2(_,xb,mu2)/xb for _ in self.iflav])
        

        # g,u,ub,d,db,s,sb,c,cb,b,bb
        # 0 1  2 3  4 5  6 7  8 9,10
        Nf=conf['alphaS'].get_Nf(mu2)

        iNf=np.zeros(11)
        iNf[1:7]=np.ones(6)
        if Nf>3: iNf[7:9]=np.ones(2)
        if Nf>4: iNf[9:11]=np.ones(2)

        lumZA  =np.sum(iNf[[1,2,3,4,5,6,7,8,9,10]]*self.eq2[[1,2,3,4,5,6,7,8,9,10]]*fA_pionZ[[0]]*fB_tungsten0[[1,2,3,4,5,6,7,8,9,10]])

        Rqg1=(za**2+(1-za)**2)/za/pT**2*lumZA/(1-xb)

        return Rqg1*prefactor

    def get_R1(self,Q,pT,s,y,chan='qqb'):
        rtau = Q/s**0.5
        xa = rtau*np.exp(y)
        xb = rtau*np.exp(-y)
        if chan=='qqb': dzadzb = lambda za,zb: self.Rqqb1zazb(Q,pT,s,y,za,zb)
        elif chan=='qg': dzadzb = lambda za,zb: self.Rqg1zazb(Q,pT,s,y,za,zb)
        elif chan=='gq': dzadzb = lambda za,zb: self.Rgq1zazb(Q,pT,s,y,za,zb)
        dza = lambda za: fixed_quad(lambda zb: np.vectorize(dzadzb)(za,zb),xb,1,n=5)[0]
        return fixed_quad(np.vectorize(dza),xa,1,n=5)[0]

    def get_asy_3diff(self,Q,pT,s,y):
        #--This computes dsig/dQ2/dy/dpT2
        qqb = self.get_R1(Q,pT,s,y,chan='qqb')
        qg  = self.get_R1(Q,pT,s,y,chan='qg')
        gq  = self.get_R1(Q,pT,s,y,chan='gq')

        return qqb+qg+gq

    def get_asy_Q_term(self,Q,pT,s,ymax):
        #--set up the integration here

        dsig=lambda y: self.get_asy_3diff(Q,pT,s,y)
        sig = fixed_quad(np.vectorize(dsig),-ymax,ymax,n=10)[0]

        jac = 2*Q * 2*pT
        return sig * jac * 0.3893793656*1e-27

    def get_asy_xF_term(self,pT,xF,s,Qmin,Qmax):
        #--set up the integration here
        mT=lambda Q: np.sqrt(pT**2+Q**2)
        y=lambda Q: np.arcsinh(xF*s**0.5/2/mT(Q))
        xFjac = lambda Q: s**0.5/mT(Q)/np.cosh(y(Q))

        dsig=lambda Q: xFjac(Q)*2*Q*self.get_asy_3diff(Q,pT,s,y(Q))
        sig = fixed_quad(np.vectorize(dsig),Qmin,Qmax,n=5)[0]

        jac = 2*pT
        return sig * jac * 0.3893793656*1e-27



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
    conf['dy_pT_scale']=0
    #pman=parman.PARMAN()
    #pman.set_qpdf_params() 
    conf['datasets']={}
    conf['datasets']['pion_qT']={}
    conf['datasets']['pion_qT']['xlsx']={}
    #conf['datasets']['pion_qT']['xlsx'][1001]='pion_qT/expdata/1001.xlsx'
    #conf['datasets']['pion_qT']['xlsx'][1002]='pion_qT/expdata/1002.xlsx'
    conf['datasets']['pion_qT']['xlsx'][10001]='pion_qT/expdata/10001.xlsx'
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

            elif tabs[idx]['obs'][0] == 'd2sigma/dpTdm (asy)':
                Q    = tabs[idx]['Q'][i]
                pT   = tabs[idx]['pT'][i]
                s    = tabs[idx]['s'][i]
                ymax = tabs[idx]['ymax'][i]
                thy = conf['pion_qT'].get_asy_Q_term(Q,pT,s,ymax)

            print idx,i,thy





