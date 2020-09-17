#!/usr/bin/env python
import os,time
import sys
import numpy as np
import pandas as pd
from scipy.integrate        import quad
from scipy.interpolate      import griddata
from obslib.idis.reader     import READER
from obslib.idis.theory     import STFUNCS
from obslib.idis.residuals  import RESIDUALS
from tools.config           import conf
from qcdlib import aux, eweak, pdf2, alphaS, mellin, ht0, ht1, ht2,offshell0, offshell1, offshell2, offshell3
from nuclib import deuterium,helium
from qcdlib.ht0 import T4
from qcdlib.ht1 import T4
from tools.tools import lprint

#--matplotlib
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('text',usetex=True)
import pylab  as py
from matplotlib.lines import Line2D

def get_grid():

    #rs    = 320.0  #--hera
    rs    = 800.0
    s     = rs**2
    Q2min = 1 #--GeV^2
    Q2max = s
    xmin  = Q2min/s
    X1    = 10**np.linspace(np.log10(xmin),-1,10)
    X2    = np.linspace(0.101,1,10)
    X     = np.append(X1,X2)
    #Q21    = np.linspace(Q2min,50,10)
    Q21    = 10**np.linspace(np.log10(Q2min),np.log10(50),5)
    Q22    = 10**np.linspace(np.log10(51),np.log10(Q2max),10)
    Q2     = np.append(Q21,Q22)
    W2max = 0.0
    M2    = conf['aux'].M2

    _X, _Q2=[],[]
    for x in X:
        for q2 in Q2:

            q2max=x*(s-M2)
            W2=M2+q2/x*(1-x)

            if q2>q2max: continue
            if W2<W2max: continue
            _X.append(x)
            _Q2.append(q2)

    return _X,_Q2

def get_data_kin():

    conf['datasets']={}
    conf['datasets']['idis']={}
    conf['datasets']['idis']['xlsx']={}
    #------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['xlsx'][10010]='idis/expdata/10010.xlsx' # proton   | F2            | SLAC                  
    conf['datasets']['idis']['xlsx'][10016]='idis/expdata/10016.xlsx' # proton   | F2            | BCDMS                 
    conf['datasets']['idis']['xlsx'][10020]='idis/expdata/10020.xlsx' # proton   | F2            | NMC                   
    conf['datasets']['idis']['xlsx'][10026]='idis/expdata/10026.xlsx' # proton   | sigma red     | HERA II NC e+ (1)     
    conf['datasets']['idis']['xlsx'][10027]='idis/expdata/10027.xlsx' # proton   | sigma red     | HERA II NC e+ (2)     
    conf['datasets']['idis']['xlsx'][10028]='idis/expdata/10028.xlsx' # proton   | sigma red     | HERA II NC e+ (3)     
    conf['datasets']['idis']['xlsx'][10029]='idis/expdata/10029.xlsx' # proton   | sigma red     | HERA II NC e+ (4)     
    conf['datasets']['idis']['xlsx'][10030]='idis/expdata/10030.xlsx' # proton   | sigma red     | HERA II NC e-         
    conf['datasets']['idis']['xlsx'][10003]='idis/expdata/10003.xlsx' # proton   | sigma red     | JLab Hall C (E00-106) 
    conf['datasets']['idis']['xlsx'][10007]='idis/expdata/10007.xlsx' # proton   | sigma red     | HERMES                
    #-------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['xlsx'][10011]='idis/expdata/10011.xlsx' # deuteron | F2            | SLAC                  
    conf['datasets']['idis']['xlsx'][10017]='idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS                 
    conf['datasets']['idis']['xlsx'][10021]='idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC                   
    conf['datasets']['idis']['xlsx'][10006]='idis/expdata/10006.xlsx' # deuteron | F2            | HERMES               
    conf['datasets']['idis']['xlsx'][10002]='idis/expdata/10002.xlsx' # deuteron | F2            | JLab Hall C (E00-106) 
    conf['datasets']['idis']['xlsx'][10033]='idis/expdata/10033.xlsx' # n/d      | F2n/F2d       | BONUS                 

    conf['datasets']['idis']['filters']=[]
    conf['datasets']['idis']['filters'].append("Q2>1.0") 
    conf['datasets']['idis']['filters'].append("W2>4.0") 

    tabs=READER().load_data_sets('idis')

    X,Q2,alpha=[],[],[]
    for idx in tabs:
        X.extend(tabs[idx]['X'])
        Q2.extend(tabs[idx]['Q2'])
        K=[_ for _ in tabs[idx].keys() if _.endswith('_u')]
        alpha2=np.zeros(len(tabs[idx]['X']))
        for _ in K: alpha2+=tabs[idx][_]**2
        alpha.extend(np.sqrt(alpha2)/tabs[idx]['value'])
    return X,Q2,np.array(alpha)

def test():

    conf['order']        = 'NLO'
    conf['Q20']          = 1.0
    conf['dglap mode']   ='truncated'
    conf['aux']          = aux.AUX()
    conf['mellin']       = mellin.MELLIN(npts=4)
    conf['dmellin']      = mellin.DMELLIN(nptsN=4,nptsM=4)
    conf['alphaS']       = alphaS.ALPHAS()
    conf['eweak']        = eweak.EWEAK()
    conf['pdf']          = pdf2.PDF()
    conf['ht']           = False
    conf['nuc']          = False
    conf['offshell']     = False
    conf['tmc']          = False
    conf['path2idistab'] = '%s/grids/grids-idis/notmc4'%os.environ['FITPACK']


    stfuncs=STFUNCS()

    X0,Q20     = get_grid()
    #GX,GLQ2    = np.meshgrid(X0,np.log(Q20))
    X,Q2,alpha = get_data_kin()

    true0 =[]
    for i in range(len(X0)):
        x,q2=X0[i],Q20[i]
        val=stfuncs.get_FXN(x,q2,stf='F2',twist=2,nucleon='proton',Nf=None,tmc=False,precalc=False,evolve=True)
        true0.append(val)
    true0=np.array(true0)

    true1=[]
    cnt=0
    for i in range(len(X)):
        cnt+=1
        lprint('%d/%d'%(cnt,len(X)))
        x,q2=X[i],Q2[i]
        val=stfuncs.get_FXN(x,q2,stf='F2',twist=2,nucleon='proton',Nf=None,tmc=False,precalc=False,evolve=True)
        true1.append(val)
    true1=np.array(true1)


    #interp0=griddata((X0,Q20),true0,(X0,Q20), fill_value=0, method='cubic')
    #interp1=griddata((X0,Q20),true0,(X,Q2), fill_value=0, method='cubic')

    interp0=griddata((np.log(X0),np.log(Q20)),true0,(np.log(X0),np.log(Q20)), fill_value=0, method='cubic')
    interp1=griddata((np.log(X0),np.log(Q20)),true0,(np.log(X) ,np.log(Q2)) , fill_value=0, method='cubic')

    alpha_=alpha*true1
    chi2=((true1-interp1)/alpha_)**2
    data=pd.DataFrame({'X':X,'Q2':Q2,'chi2':chi2})

    dA=data.query('chi2<0.1')
    dB=data.query('chi2>0.1')
    dC=data.query('chi2>1')

    nrows,ncols=3,2
    fig = py.figure(figsize=(ncols*5,nrows*3))

    #-----------------
    ax=py.subplot(nrows,ncols,1)
    ax.plot(X,Q2,'b.')
    ax.plot(X0,Q20,'r.')
    ax.semilogx()
    ax.semilogy()
    ax.text(0.1,0.8,r'${\rm grid~points}=%d$'%len(X0),transform=ax.transAxes,size=20)
    ax.text(0.1,0.7,r'${\rm exp.~points}=%d$'%len(X),transform=ax.transAxes,size=20)
    ax.set_ylabel(r'$Q^2$',size=20)

    #-----------------
    ax=py.subplot(nrows,ncols,2)
    ax.plot(X,Q2,'b.')
    ax.plot(X0,Q20,'r.')
    ax.semilogy()

    #-----------------
    ax=py.subplot(nrows,ncols,3)
    ax.plot(X0,Q20,'r.')
    ax.plot(dB.X,dB.Q2,'m.',label=r'$\chi^2>0.1$')
    ax.plot(dC.X,dC.Q2,'b.',label=r'$\chi^2>1$')
    ax.legend()
    ax.semilogx()
    ax.semilogy()
    ax.set_ylabel(r'$Q^2$',size=20)

    #-----------------
    ax=py.subplot(nrows,ncols,4)
    ax.plot(X0,Q20,'r.')
    ax.plot(dB.X,dB.Q2,'m.',label=r'$\chi^2>0.1$')
    ax.plot(dC.X,dC.Q2,'b.',label=r'$\chi^2>1$')
    #ax.semilogx()
    ax.semilogy()

    #-----------------
    ax=py.subplot(nrows,ncols,5)
    ax.plot(dA.X,dA['chi2'],'y.',label=r'$\chi^2<0.1$')
    ax.plot(dB.X,dB['chi2'],'m.',label=r'$\chi^2>0.1$')
    ax.plot(dC.X,dC['chi2'],'b.',label=r'$\chi^2>1$')
    ax.axhline(0.1)
    ax.axhline(1.0)
    ax.semilogx()
    ax.semilogy()
    ax.set_ylim(0.001,None)
    ax.set_ylabel(r'$\chi^2$',size=20)
    ax.set_xlabel(r'$x_{\rm bj}$',size=20)

    #-----------------
    ax=py.subplot(nrows,ncols,6)
    ax.plot(dA.X,dA['chi2'],'y.',label=r'$\chi^2<0.1$')
    ax.plot(dB.X,dB['chi2'],'m.',label=r'$\chi^2>0.1$')
    ax.plot(dC.X,dC['chi2'],'b.',label=r'$\chi^2>1$')
    ax.axhline(0.1)
    ax.axhline(1.0)
    ax.semilogy()
    ax.set_ylim(0.001,None)
    ax.set_xlabel(r'$x_{\rm bj}$',size=20)
    #-----------------
    py.tight_layout()
    py.savefig('test.pdf')

    print('\nSummary')
    print('grid points=%d'%len(X0))
    print('exp  points=%d'%len(X))

    return 


if __name__=="__main__":

    test()




















