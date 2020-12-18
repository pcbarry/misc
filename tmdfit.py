#!/usr/bin/env python
import os,sys
import numpy as np
import pandas as pd
from scipy.integrate import quad,fixed_quad
from scipy.optimize import leastsq
from scipy.special import jv

from tools.config import conf,load_config
from tools.tools import load,save,checkdir,lprint

from timeit import default_timer as timer

from obslib.pion_qT.reader import READER
from qcdlib import aux,mellin,tmdpdf

import copy
import lhapdf


def get_collinear(hadA,hadB,wdir):
	load_config('%s/input.py'%wdir)
	conf['pdf-pion']=lhapdf.mkPDF(hadA)
	conf['pdf-W']=lhapdf.mkPDF(hadB)

def get_Wtilde(bT,Q,xA,xB,par,gjoption=1):

    zeta=Q**2
    mub=tmdpi.get_mub(bT)
    mu=Q

    if (xA,xB,mub) not in conf['A storage']: 
        A=0
        flavor=[['d','db'],['u','ub'],['s','sb']]
        if conf['aux'].mc2<=mub: flavor.append(['c','cb'])
        if conf['aux'].mb2<=mub: flavor.append(['b','bb'])
        cnt=1
        for flav in flavor:
            if flav[0]=='u' or flav[0]=='c': eq2=4.0/9.0
            else: eq2=1.0/9.0
            fA=conf['pdf-pion'].xfxQ2(cnt,xA,mub**2)/xA
            fB=conf['pdf-W'].xfxQ2(-cnt,xB,mub**2)/xB
            A+=eq2*fA*fB
            fA=conf['pdf-pion'].xfxQ2(-cnt,xA,mub**2)/xA
            fB=conf['pdf-W'].xfxQ2(cnt,xB,mub**2)/xB
            A+=eq2*fA*fB
            cnt+=1
        conf['A storage'][(xA,xB,mub)]=A

    A=conf['A storage'][(xA,xB,mub)]

    B=np.exp(2*tmdpi.get_LB(bT,mu,zeta))

    gKpiece=(Q**2/conf['aux'].mc2)**(-par[0]*(1-np.exp(-tmdpi.CF*tmdpi.get_alpha(mub**2)*bT**2/(np.pi*par[0]*tmdpi.bTmax**2))))

    if gjoption==1:
        C1A=np.exp(-bT**2*par[1])* par[2]*xA**par[3]*(1-xA)**par[4]
        C1B=np.exp(-bT**2*par[1])* par[2]*xB**par[3]*(1-xB)**par[4]

    if gjoption==2:
        C1A=np.exp(-bT**2*par[1])* (par[2]*xA)**(-par[1]*par[3]*bT**2)
        C1B=np.exp(-bT**2*par[1])* (par[2]*xB)**(-par[1]*par[3]*bT**2)

    if gjoption==3:
        C1A=np.exp(-bT**2*par[1])* par[2]*xA**(-par[1]*par[3]*bT**2)
        C1B=np.exp(-bT**2*par[1])* par[2]*xB**(-par[1]*par[3]*bT**2)

    if gjoption==4:
        C1A=par[1]*np.exp(-(par[2]+(1-xA)**2*par[3])*bT**2/(1+par[4]*bT**2)**0.5)
        C1B=par[1]*np.exp(-(par[2]+(1-xB)**2*par[3])*bT**2/(1+par[4]*bT**2)**0.5)

    lam=1
    p=6
    R= np.exp(-(bT/lam)**p)

    return (A*B)*R + (1-R)*C1A*C1B*gKpiece**2

def model(npopt,par=None):

    data={'xF':[],'pT':[],'exp':[],'alpha':[],'thy':[]}

    print('params:')
    print(par)

    #--get a bin from table
    d=tab[1002].query('pT<1.0 and xF<0.8')

    N=len(d['s'].values)
    for i in range(N):
        s=d['s'].values[i]
        xF=d['xF'].values[i]
        Qmin=d['Qmin'].values[i]
        Qmax=d['Qmax'].values[i]
        pT=d['pT'].values[i]

        #--get area under Wtilde
        hc = 0.3893793721 #--GeV^2 mbarn
        hc = hc *  1e-27 #--GeV^2 cm^2
        exp=d['value'].values[i]/hc  #/(2*pT)   #--1/GeV^3
        alpha=d['stat_u'].values[i]/hc

        def integrand(Q):

            mT=(Q**2+pT**2)**0.5
            y=np.arcsinh(xF*np.sqrt(s)/2./mT)

            #--compute W in qT space
            #bTmin=tmdpi.bTgrid[0]
            #bTmax=tmdpi.bTgrid[-1]   
            bTmin=1e-4
            bTmax=1e2

            xA=(Q/np.sqrt(s))*np.exp(+y)
            xB=(Q/np.sqrt(s))*np.exp(-y)
            if xA>1: return 0
            if xB>1: return 0

            W=fixed_quad(np.vectorize(lambda bT: bT*jv(0,pT*bT)/(2*np.pi)*get_Wtilde(bT,Q,xA,xB,par,gjoption=npopt)),bTmin,bTmax,n=100)[0]

            alfa2=1/137.0**2
            born=4*np.pi**2*alfa2/9/Q**2/s
            jac =np.sqrt(s)/(2*mT*np.cosh(y))*2*pT*2*Q

            return born*jac*W

        thy=fixed_quad(np.vectorize(integrand),Qmin,Qmax,n=20)[0]
        data['xF'].append(xF)
        data['pT'].append(pT)
        data['exp'].append(exp)
        data['thy'].append(thy)
        data['alpha'].append(alpha)

    res=(np.array(data['thy'])-np.array(data['exp']))/np.array(data['alpha'])
    chi2=np.sum(res**2)
    print 'chi2: %.4f'%(chi2/N)
    return res


if __name__=='__main__':

    wdir=sys.argv[1]
    option=int(sys.argv[2])

    conf['aux']=aux.AUX()
    conf['mellin']=mellin.MELLIN()
    conf['datasets']={}
    conf['datasets']['pion_qT']={}
    conf['datasets']['pion_qT']['xlsx']={}
    #conf['datasets']['pion_qT']['xlsx'][1001]='pion_qT/expdata/1001.xlsx'
    conf['datasets']['pion_qT']['xlsx'][1002]='pion_qT/expdata/1002.xlsx'

    tab=READER().load_data_sets('pion_qT')
    for _ in tab: tab[_]=pd.DataFrame(tab[_])


    pi='JAM20PionPDFnlo'
    W='EPPS16nlo_CT14nlo_W184'
    get_collinear(pi,W,wdir)

    tmdpi=tmdpdf.QCDTMD()
    tmdW=tmdpdf.QCDTMD()

    #--storages
    conf['A storage']={}

    mod=lambda par: model(option,par)

    if option==1: PAR=[0.01,0.1425,3098.3045,6.7446,0.5593]
    elif option==2: PAR=[0.68,0.21,10,-0.6]
    elif option==3: PAR=[0.68,0.21,10,-0.6]
    elif option==4: PAR=[0.022,1,0.17,0.48,2.15]


    fit=leastsq(mod,PAR,full_output=1)
    sol=fit[0]

    res=mod(sol)

    data={}
    data['initial params']=PAR
    data['params']=sol
    data['residuals']=res
    data['chi2']=np.sum(res**2)

    checkdir('%s/fitdata'%wdir)
    save(data,'%s/fitdata/opt%i.dat'%(wdir,option))



