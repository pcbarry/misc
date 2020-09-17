#!/usr/bin/env python
import numpy as np
from tools.config       import conf
from tools.residuals    import _RESIDUALS
from obslib.dy.reader   import READER
from obslib.dy.theory   import DY
from obslib.dy.dymellin import DYMELLIN
from obslib.dy.piontheory import DY_PION
from obslib.dy.piontheoryplus import DY_PION_PLUS
from obslib.dy.piontheorynlo import DY_PION_NLO

class RESIDUALS(_RESIDUALS):

  def __init__(self): 
    self.reaction='dy'
    self.tabs=conf['dy tabs']
    self.setup()

  def _get_theory(self,entry):
    k,i=entry
    Q2    =self.tabs[k]['Q2'][i]
    Jac   =self.tabs[k]['Jac'][i]
    units =self.tabs[k]['Units'][i]
    reaction =self.tabs[k]['reaction'][i]
    thy=conf['dy'].get_mxsec(k,i,Q2,reaction) * Jac * units
    return thy

  def gen_report(self,verb=1,level=1):
    """
    verb = 0: Do not print on screen. Only return list of strings
    verv = 1: print on screen the report
    level= 0: only the total chi2s
    level= 1: include point by point 
    """
    L=[]

    L.append('reaction: dy')
    msg ='%7s'%'idx'
    msg+='%10s'%'tar'
    msg+='%10s'%'col'
    msg+='%5s'%'npts'
    msg+='%10s'%'chi2'
    msg+='%10s'%'chi2-npts'
    msg+='%10s'%'rchi2'
    msg+='%10s'%'nchi2'
    L.append(msg)
    for idx in self.tabs:
      res =self._get_residuals(idx)
      rres=self._get_rres(idx)
      nres=self._get_nres(idx)
      
      chi2=np.sum(res**2)
      rchi2=np.sum(rres**2)
      nchi2=nres**2
      reaction=self.tabs[idx]['reaction'][0]
      col=self.tabs[idx]['col'][0].split()[0]
      npts=res.size
      msg ='%7d'%idx
      msg+='%10s'%reaction
      msg+='%10s'%col
      msg+='%5d'%npts
      msg+='%10.2f'%chi2
      msg+='%10.2f'%(chi2-npts)
      msg+='%10.2f'%rchi2
      msg+='%10.2f'%nchi2
      L.append(msg)

    if level==1:
      L.append('-'*100)  

      msg ='col=%7s,  '
      msg+='tar=%7s,  '
      msg+='obs=%7s,  '
      msg+='xF=%10.3e,  '
      msg+='Y=%10.3e,  '
      msg+='Q2=%10.3e,  '
      msg+='exp=%10.3e,  ' 
      msg+='alpha=%10.3e,  ' 
      msg+='thy=%10.3e,  ' 
      msg+='shift=%10.3e,  ' 
      msg+='chi2=%10.3f,  '
      msg+='norm=%10.3f,  '

      for k in self.tabs:
        if k in conf['datasets']['dy']['norm']:
          norm=conf['datasets']['dy']['norm'][k]['value']
        else:
          norm=1
        for i in range(len(self.tabs[k]['value'])):
          xF    = self.tabs[k]['xF'][i]
          Y     = self.tabs[k]['Y'][i]
          obs   = self.tabs[k]['obs'][i]
          Q2    = self.tabs[k]['Q2'][i]
          res   = self.tabs[k]['residuals'][i]
          thy   = self.tabs[k]['thy'][i]
          exp   = self.tabs[k]['value'][i]
          alpha = self.tabs[k]['alpha'][i]
          rres  = self.tabs[k]['r-residuals'][i]
          col   = self.tabs[k]['col'][i]
          shift = self.tabs[k]['shift'][i]
          reaction = self.tabs[k]['reaction'][i]
          if res<0: chi2=-res**2
          else: chi2=res**2
          L.append(msg%(col,reaction,obs,xF,Y,Q2,exp,alpha,thy,shift,chi2,norm))

    if verb==0:
      return L
    elif verb==1:
      for l in L: print l
      return L

class RESIDUALS_PROTONS(_RESIDUALS):

  def __init__(self):
    self.reaction='dy-pion'
    self.tabs=conf['dy-pion tabs']
    self.setup()

  def _get_theory(self,entry): 
    k,i=entry
    Q2    =self.tabs[k]['Q2'][i]
    Jac   =self.tabs[k]['Jac'][i]
    units =self.tabs[k]['Units'][i]
    obs   =self.tabs[k]['obs'][i]
    if obs == 'dsig/drtau/dxF (valence)':
        conf['pdf'].evolve(Q2)
        for f in ['g','ub','db','s','sb','c','cb','b','bb']:
            conf['pdf'].storage[Q2][f]=np.zeros(len(conf['pdf'].storage[Q2][f]))
        thy1=conf['dy-pion'].get_xsec_mell_hybrid(k,i,Q2) * Jac * units
        del conf['pdf'].storage[Q2]

        conf['pdf'].evolve(Q2)
        for f in ['g','u','d','s','sb','c','cb','b','bb']:
            conf['pdf'].storage[Q2][f]=np.zeros(len(conf['pdf'].storage[Q2][f]))
        thy2=conf['dy-pion'].get_xsec_mell_hybrid(k,i,Q2) * Jac * units
        del conf['pdf'].storage[Q2]
        thy = thy1 - thy2

    elif obs == 'dsig/drtau/dxF (sea)':
        conf['pdf'].evolve(Q2)
        conf['pdf'].storage[Q2]['g']=np.zeros(len(conf['pdf'].storage[Q2]['g']))
        conf['pdf'].storage[Q2]['u']=conf['pdf'].storage[Q2]['ub']
        conf['pdf'].storage[Q2]['d']=conf['pdf'].storage[Q2]['db']
        thy=conf['dy-pion'].get_xsec_mell_hybrid(k,i,Q2) * Jac * units
        del conf['pdf'].storage[Q2]

    elif obs == 'dsig/drtau/dxF (gluon)':
        conf['pdf'].evolve(Q2)
        for f in ['u','ub','d','db','s','sb','c','cb','b','bb']:
            conf['pdf'].storage[Q2][f]=np.zeros(len(conf['pdf'].storage[Q2][f]))
        thy=conf['dy-pion'].get_xsec_mell_hybrid(k,i,Q2) * Jac * units
        del conf['pdf'].storage[Q2]

    #elif obs == 'dsig/drtau/dxF (plus)':
    #    thy = conf['dy-pion plus'].get_xsec_mell_hybrid(k,i,Q2) * Jac * units

    #elif obs == 'dsig/drtau/dxF (NLO)':
    #    thy = conf['dy-pion nlo'].get_xsec_mell_hybrid(k,i,Q2) * Jac * units

    #elif obs == 'dsig/drtau/dxF (LO)':
    #    #dy=DY_PION_NLO()
    #    #dy.load_melltab_hybrid()
    #    conf['dy-pion lo'].iord=0
    #    thy = conf['dy-pion lo'].get_xsec_mell_hybrid(k,i,Q2) * Jac * units

    else:
        thy=conf['dy-pion'].get_xsec_mell_hybrid(k,i,Q2) * Jac * units
    #if type(conf['pdf-pion parametrization'])==str and 'resum' in conf['pdf-pion parametrization']:
    #  S=self.tabs[k]['S'][i]
    #  Y=self.tabs[k]['Y'][i]
    #  thy+=conf['resummed dy-pion'].get_interpolated(S,Y,Q2,Q2) * Jac * units
    return thy

  def gen_report(self,verb=1,level=1):
    """
    verb = 0: Do not print on screen. Only return list of strings
    verv = 1: print on screen the report
    level= 0: only the total chi2s
    level= 1: include point by point 
    """
    L=[]

    L.append('reaction: dy')
    msg ='%7s'%'idx'
    #msg+='%10s'%'reaction'
    msg+='%10s'%'col'
    msg+='%5s'%'npts'
    msg+='%10s'%'chi2'
    msg+='%10s'%'chi2-npts'
    msg+='%10s'%'rchi2'
    msg+='%10s'%'nchi2'
    L.append(msg)
    for idx in self.tabs:
      res =self._get_residuals(idx)
      rres=self._get_rres(idx)
      nres=self._get_nres(idx)
      
      chi2=np.sum(res**2)
      rchi2=np.sum(rres**2)
      nchi2=nres**2
      #reaction=self.tabs[idx]['reaction'][0]
      col=self.tabs[idx]['col'][0].split()[0]
      npts=res.size
      msg ='%7d'%idx
      #msg+='%10s'%reaction
      msg+='%10s'%col
      msg+='%5d'%npts
      msg+='%10.2f'%chi2
      msg+='%10.2f'%(chi2-npts)
      msg+='%10.2f'%rchi2
      msg+='%10.2f'%nchi2
      L.append(msg)

    if level==1:
      L.append('-'*100)  

      msg ='col=%7s,  '
      msg+='tar=%7s,  '
      msg+='obs=%7s,  '
      msg+='xF=%10.3e,  '
      msg+='Y=%10.3e,  '
      msg+='Q2=%10.3e,  '
      msg+='exp=%10.3e,  ' 
      msg+='alpha=%10.3e,  ' 
      msg+='thy=%10.3e,  ' 
      msg+='shift=%10.3e,  ' 
      msg+='chi2=%10.3f,  '
      msg+='norm=%10.3f,  '

      for k in self.tabs:
        if k in conf['datasets']['dy']['norm']:
          norm=conf['datasets']['dy']['norm'][k]['value']
        else:
          norm=1
        for i in range(len(self.tabs[k]['value'])):
          xF    = self.tabs[k]['xF'][i]
          Y     = self.tabs[k]['Y'][i]
          obs   = self.tabs[k]['obs'][i]
          Q2    = self.tabs[k]['Q2'][i]
          res   = self.tabs[k]['residuals'][i]
          thy   = self.tabs[k]['thy'][i]
          exp   = self.tabs[k]['value'][i]
          alpha = self.tabs[k]['alpha'][i]
          rres  = self.tabs[k]['r-residuals'][i]
          col   = self.tabs[k]['col'][i]
          shift = self.tabs[k]['shift'][i]
          #reaction = self.tabs[k]['reaction'][i]
          if res<0: chi2=-res**2
          else: chi2=res**2
          #L.append(msg%(col,reaction,obs,xF,Y,Q2,exp,alpha,thy,shift,chi2,norm))
          L.append(msg%(col,obs,xF,Y,Q2,exp,alpha,thy,shift,chi2,norm))

    if verb==0:
      return L
    elif verb==1:
      for l in L: print l
      return L
    """
    verb = 0: Do not print on screen. Only return list of strings
    verb = 1: print on screen the report
    level= 0: only the total chi2s
    level= 1: include point by point
    """
    L=[]

    L.append('reaction: dy-pion')

    L.append('%7s %10s %10s %5s %10s %10s %10s %10s %10s'%('idx','tar','col','npts','chi2','chi2-npts','chi2/npts','rchi2','nchi2'))
    for k in self.tabs:
      res =self.tabs[k]['residuals']
      rres=self._get_rres(k)

if __name__=='__main__':
  
    import os
    from qcdlib import alphaS,aux,eweak,mellin,pdf0
    from qcdlib import aux,mellin,alphaS,eweak
    from reader import READER  
    import time

    conf['Q20']   = 1.0
    conf['alphaSmode']='backward'
    conf['dglap mode']='truncated'
    conf['order']='NLO'
    conf['scheme']='ZMVFS'
    conf['path2dytab']='%s/grids/dytab'%os.environ['FITPACK']
    conf['aux']=aux.AUX()
    conf['mellin']=mellin.MELLIN(npts=4)
    conf['alphaS']=alphaS.ALPHAS()
    conf['eweak']=eweak.EWEAK()
    conf['pdfA']=pdf0.PDF()
    conf['pdfB']=conf['pdfA']
    conf['datasets']={}
    conf['datasets']['dy']={}
    conf['datasets']['dy']['xlsx']={}
    conf['datasets']['dy']['xlsx'][10001]='dy/expdata/10001.xlsx'
    conf['datasets']['dy']['xlsx'][10002]='dy/expdata/10002.xlsx'
    conf['datasets']['dy']['norm']={}
    conf['datasets']['dy']['filters']=[]
    conf['datasets']['dy']['filters'].append("Q2>1") 
    conf['dy tabs']=READER().load_data_sets('dy')
    conf['dy']=DYMELLIN()
    conf['dy'].load_melltab()

    residuals=RESIDUALS()
    t1=time.time()
    residuals.get_residuals()
    print(residuals.get_residuals())
    t2=time.time()
    print t2-t1

