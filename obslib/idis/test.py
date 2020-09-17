#!/usr/bin/env python
import os,time
import sys
import numpy as np
from scipy.integrate import quad
from obslib.idis.reader     import READER
from obslib.idis.theory     import STFUNCS
#from obslib.idis.theory_ML     import STFUNCS as STFUNCS_ML
from obslib.idis.residuals  import RESIDUALS
from tools.config           import conf
from qcdlib import aux, eweak, pdf2, alphaS, mellin, ht0, ht1, ht2,offshell0, offshell1, offshell2, offshell3
from nuclib import deuterium,helium
from qcdlib.ht0 import T4
from qcdlib.ht1 import T4

def test_theory():

    conf['order'] = 'NLO'
    conf['Q20']   = 1.0
    conf['dglap mode']='truncated'
    conf['aux']   = aux.AUX()
    conf['mellin']= mellin.MELLIN(npts=4)
    conf['dmellin']=mellin.DMELLIN(nptsN=4,nptsM=4)
    conf['alphaS']= alphaS.ALPHAS()
    conf['eweak'] = eweak.EWEAK()
    conf['pdf']   = pdf2.PDF()
    conf['ht']=True
    conf['ht parametrization']=2
    if conf['ht parametrization']==0:
      conf['t4F2']  = ht0.T4()
      conf['t4FL']  = ht0.T4()
      conf['t4F3']  = ht0.T4()
    if conf['ht parametrization']==1:
      conf['t4F2']  = ht1.T4()
      conf['t4FL']  = ht1.T4()
      conf['t4F3']  = ht1.T4()
    if conf['ht parametrization']==2:
      conf['t4F2']  = ht2.T4()
      conf['t4FL']  = ht2.T4()
      conf['t4F3']  = ht2.T4()
    conf['offshell parametrization']=2
    if conf['offshell parametrization']==0: conf['F2off']  = offshell0.OFFSHELL()
    if conf['offshell parametrization']==1: conf['F2off']  = offshell1.OFFSHELL()
    if conf['offshell parametrization']==2: conf['F2off']  = offshell2.OFFSHELL()
    if conf['offshell parametrization']==3: conf['F2off']  = offshell3.OFFSHELL()
    conf['nuc']    = True
    conf['offshell'] = False
    conf['tmc'] = False
    if   conf['tmc']==False: conf['path2idistab'] = '%s/grids/grids-idis/notmc4'%os.environ['FITPACK']
    elif conf['tmc']=='GP' : conf['path2idistab'] = '%s/grids/grids-idis/gp4'%os.environ['FITPACK']
    elif conf['tmc']=='AOT': conf['path2idistab'] = '%s/grids/grids-idis/aot4'%os.environ['FITPACK']

    stfuncs=STFUNCS()
    x=0.083
    Q2=2.39
    #print stfuncs.get_FXN(x,Q2,stf='F2',twist=4,nucleon='proton',Nf=None,tmc=False,precalc=False,evolve=True)
    #print stfuncs.get_F2D(x,Q2,twist=4,Nf=None,tmc=False,nuc=False,offshell=False,evolve=True)
    print stfuncs.get_F2D(x,Q2,twist=4,Nf=None,tmc=False,nuc=False,offshell=False,evolve=False)
    #print stfuncs.get_FLD(x,Q2,twist=4,Nf=None,tmc=False,nuc=True,offshell=False,evolve=True)


def test_theory_ML():

    conf['order'] = 'NLO'
    conf['Q20']   = 1.0
    conf['dglap mode']='truncated'
    conf['aux']   = aux.AUX()
    conf['mellin']= mellin.MELLIN(npts=4)
    conf['dmellin']=mellin.DMELLIN(nptsN=4,nptsM=4)
    conf['alphaS']= alphaS.ALPHAS()
    conf['eweak'] = eweak.EWEAK()
    conf['pdf']   = pdf2.PDF()
    conf['ht']=True
    conf['ht parametrization']=2
    if conf['ht parametrization']==0:
      conf['t4F2']  = ht0.T4()
      conf['t4FL']  = ht0.T4()
      conf['t4F3']  = ht0.T4()
    if conf['ht parametrization']==1:
      conf['t4F2']  = ht1.T4()
      conf['t4FL']  = ht1.T4()
      conf['t4F3']  = ht1.T4()
    conf['F2off']  = ht0.T4()
    conf['FLoff']  = ht0.T4() 
    conf['F3off']  = ht0.T4()
    conf['nuc']    = True
    conf['offshell'] = False
    conf['tmc'] = False
    if   conf['tmc']==False: conf['path2idistab'] = '/work/JAM/ccocuzza/ml4jam/idis/models/grids-idis/notmc4'
    elif conf['tmc']=='GP' : conf['path2idistab'] = '%s/grids/grids-idis/gp4'%os.environ['FITPACK']
    elif conf['tmc']=='AOT': conf['path2idistab'] = '%s/grids/grids-idis/aot4'%os.environ['FITPACK']

    stfuncs=STFUNCS()
    x=0.083
    Q2=2.39
    #print stfuncs.get_FXN(x,Q2,stf='F2',twist=4,nucleon='proton',Nf=None,tmc=False,precalc=False,evolve=True)
    #print stfuncs.get_F2D(x,Q2,twist=4,Nf=None,tmc=False,nuc=False,offshell=False,evolve=True)
    print stfuncs.get_F2D(x,Q2,twist=4,Nf=None,tmc=False,nuc=True,offshell=False,evolve=True)
    #print stfuncs.get_FLD(x,Q2,twist=4,Nf=None,tmc=False,nuc=True,offshell=False,evolve=True)
    #print hqstfuncs.get_FXHQ(x,Q2,stf='F2',hq='c',twist=4,nucleon='proton',Nf=None,tmc=False,evolve=True)
#    print hqstfuncs.get_FXHQ(x,Q2,stf='F2',hq='c',twist=2,nucleon='proton',Nf=None,tmc=False,evolve=True)


def test_offshell():
    #--calculate offshell through integral and compare to moments method
    conf['order'] = 'NLO'
    conf['Q20']   = 1.3**2
    conf['dglap mode']='truncated'
    conf['aux']   = aux.AUX()
    conf['mellin']= mellin.MELLIN(npts=4)
    conf['dmellin']=mellin.DMELLIN(nptsN=4,nptsM=4)
    conf['alphaS']= alphaS.ALPHAS()
    conf['eweak'] = eweak.EWEAK()
    conf['pdf']   = pdf2.PDF()
    conf['dsmf']  = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])
    conf['ht']=False
    conf['offshell parametrization']=2
    if conf['offshell parametrization']==0: conf['F2off']  = offshell0.OFFSHELL()
    if conf['offshell parametrization']==1: conf['F2off']  = offshell1.OFFSHELL()
    if conf['offshell parametrization']==2: conf['F2off']  = offshell2.OFFSHELL()
    if conf['offshell parametrization']==3: conf['F2off']  = offshell3.OFFSHELL()
    conf['nuc']    = True
    conf['offshell'] = True
    conf['tmc'] = False
    if   conf['tmc']==False: conf['path2idistab'] = '%s/grids/grids-idis/notmc4'%os.environ['FITPACK']
    elif conf['tmc']=='GP' : conf['path2idistab'] = '%s/grids/grids-idis/gp4'%os.environ['FITPACK']
    elif conf['tmc']=='AOT': conf['path2idistab'] = '%s/grids/grids-idis/aot4'%os.environ['FITPACK']

    stfuncs=STFUNCS()


    x = 0.8
    Q2 = 13.88

    ymax = 1.997635579 #mD/mN

    #--onshell part
    smear = lambda y: conf['dsmf'].get_fXX('f22','onshell',x,Q2,y)
    stf = lambda y: stfuncs.get_FXN(x/y,Q2,twist=2,nucleon='proton',tmc=conf['tmc'])\
                   +stfuncs.get_FXN(x/y,Q2,twist=2,nucleon='neutron',tmc=conf['tmc'])

    integrand = lambda y: smear(y)*stf(y)

    integral_result = quad(integrand,x,ymax,full_output=1)[0]

    print('Integral result (onshell): %6.5f'%integral_result)
    print('Mellin result (onshell):   %6.5f'%stfuncs.get_F2D(x,Q2,twist=2,nuc=True,tmc=conf['tmc'],offshell=False))

    #--offshell part
    integral_result=0
    smear = lambda y: conf['dsmf'].get_fXX('f22','offshell',x,Q2,y)
    for nucleon in ['proton','neutron']:
      if conf['offshell parametrization']==0 or conf['offshell parametrization']==1:
        norm, a, b, c, d = conf['F2off'].params[nucleon]
        off = lambda y: norm*((x/y)**a)*((1-(x/y))**b)*(1+c*((x/y)**0.5)+d*(x/y))
      if conf['offshell parametrization']==2 or conf['offshell parametrization']==3:
        norm, x0, x1 = conf['F2off'].params[nucleon]
        off = lambda y: norm*((x/y) - x0)*((x/y) - x1)*(1 + x0 - (x/y))
      stf = lambda y: stfuncs.get_FXN(x/y,Q2,stf='F2',twist=2,nucleon=nucleon,tmc=conf['tmc'])
      integrand = lambda y: smear(y) * off(y) * stf(y)
      integral_result += quad(integrand,x,ymax,full_output=1)[0]


    print('Integral result (offshell): %6.5f'%integral_result)
    print('Mellin result (offshell):   %6.5f'%stfuncs.get_F2D(x,Q2,twist=2,nuc=True,tmc=conf['tmc'],offshell=True))

def test_reader():

    conf['datasets']={}

    conf['datasets']['idis']={}
    conf['datasets']['idis']['xlsx']={}
    conf['datasets']['idis']['xlsx'][10010]='idis/expdata/10010.xlsx' # proton   | F2            | SLAC                  
    conf['datasets']['idis']['xlsx'][10016]='idis/expdata/10016.xlsx' # proton   | F2            | BCDMS                 
    conf['datasets']['idis']['xlsx'][10020]='idis/expdata/10020.xlsx' # proton   | F2            | NMC                   
    conf['datasets']['idis']['xlsx'][10026]='idis/expdata/10026.xlsx' # proton   | sigma red     | HERA II NC e+ (1)     
    conf['datasets']['idis']['xlsx'][10027]='idis/expdata/10027.xlsx' # proton   | sigma red     | HERA II NC e+ (2)     
    conf['datasets']['idis']['xlsx'][10028]='idis/expdata/10028.xlsx' # proton   | sigma red     | HERA II NC e+ (3)     
    conf['datasets']['idis']['xlsx'][10029]='idis/expdata/10029.xlsx' # proton   | sigma red     | HERA II NC e+ (4)     
    conf['datasets']['idis']['xlsx'][10030]='idis/expdata/10030.xlsx' # proton   | sigma red     | HERA II NC e-         
    conf['datasets']['idis']['xlsx'][10003]='idis/expdata/10003.xlsx' # proton   | F2            | JLab Hall C (E00-106) 
    conf['datasets']['idis']['xlsx'][10007]='idis/expdata/10007.xlsx' # proton   | sigma red     | HERMES                

    conf['datasets']['idis']['xlsx'][10011]='idis/expdata/10011.xlsx' # deuteron | F2            | SLAC                  
    conf['datasets']['idis']['xlsx'][10017]='idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS                 
    conf['datasets']['idis']['xlsx'][10021]='idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC                   
    conf['datasets']['idis']['xlsx'][10006]='idis/expdata/10006.xlsx' # deuteron | F2            | HERMES               
    conf['datasets']['idis']['xlsx'][10002]='idis/expdata/10002.xlsx' # deuteron | F2            | JLab Hall C (E00-106) 

    conf['datasets']['idis']['xlsx'][10033]='idis/expdata/10033.xlsx' # n/d      | F2n/F2d       | BONUS                 

    conf['datasets']['idis']['filters']=[]
    conf['datasets']['idis']['filters'].append("Q2>1.0") 
    conf['datasets']['idis']['filters'].append("W2>4.0") 


    conf['aux']=aux.AUX()
    TAB=READER().load_data_sets('idis')

def test_residuals(nuc=False,tmc=False,ht=False,offshell=False,hq=False):

    print("="*50)
    print('Test setup:')
    print('nuc=',nuc)
    print('tmc=',tmc)
    print('ht=',ht)
    print('offshell=',offshell)
    print('hq=',hq)

    conf['order'] = 'NLO'
    conf['Q20']   = 1.0
    conf['dglap mode']='truncated'
    conf['aux']   = aux.AUX()
    conf['mellin']= mellin.MELLIN(npts=4)
    conf['dmellin']=mellin.DMELLIN(nptsN=4,nptsM=4)
    conf['alphaS']= alphaS.ALPHAS()
    conf['eweak'] = eweak.EWEAK()
    conf['pdf']   = pdf2.PDF()
    conf['tmc']   = tmc
    conf['ht']    = ht
    conf['hq']    = hq
    conf['nuc']   = nuc

    conf['t4F2']  = T4()
    conf['t4FL']  = T4()
    conf['t4F3']  = T4()
    conf['t4W2']  = T4()
    conf['t4WL']  = T4()
    conf['t4W3']  = T4()
  
    conf['dsmf']  = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])

    if tmc==False:
        conf['path2idistab'] = '%s/grids/grids-idis/notmc4'%os.environ['FITPACK']
    if tmc=='GP':
        conf['path2idistab'] = '%s/grids/grids-idis/gp4'%os.environ['FITPACK']
    if tmc=='AOT':
        conf['path2idistab'] = '%s/grids/grids-idis/aot4'%os.environ['FITPACK']
    conf['idis stfuncs']=STFUNCS()
  
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
    #-------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['norm']={}
    #-------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['norm'][10010]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10016]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10020]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10003]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10007]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    #-------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['norm'][10011]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10017]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10006]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10002]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10033]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}


    conf['datasets']['idis']['filters']=[]
    conf['datasets']['idis']['filters'].append("Q2>1.0") 
    conf['datasets']['idis']['filters'].append("W2>4.0") 
    conf['idis tabs']=READER().load_data_sets('idis')

    residuals=RESIDUALS()
    res,rres,nres=residuals.get_residuals()

if __name__=='__main__':

  
    test_theory()
    #test_theory_ML()

    #test_offshell()

    #test_reader()

    #test_residuals(nuc=False,tmc=False,ht=False,hq=False)
    #test_residuals(nuc=True ,tmc=False,ht=False,hq=False)

    #test_residuals(nuc=False,tmc='GP' ,ht=False,hq=False)
    #test_residuals(nuc=False,tmc='AOT',ht=False,hq=False)

    #test_residuals(nuc=True,tmc='GP' ,ht=False,hq=False)
    #test_residuals(nuc=True,tmc='AOT',ht=False,hq=False)

    #test_residuals(nuc=True,tmc='GP' ,ht=True,hq=False)
    #test_residuals(nuc=True,tmc='AOT',ht=True,hq=False)

    #test_residuals(nuc=True,tmc=False,ht=False,offshell=True,hq=False)
    #test_residuals(nuc=True,tmc=False,ht=True,offshell=True,hq=False)
